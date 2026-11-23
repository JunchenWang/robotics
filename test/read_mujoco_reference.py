"""Export URDF chains and runtime reference values for test_read_mujoco."""

import argparse
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path

import mujoco
import numpy as np


FRAME_FIXTURE = """<mujoco model="joint_frames">
  <option gravity="1.2 -2.3 -9.4"/>
  <default>
    <default class="link_geom"><geom density="600" contype="0" conaffinity="0"/></default>
  </default>
  <worldbody>
    <frame pos="0.5 0 0.1" euler="12 -15 23">
      <body name="mount" pos="0.1 0.2 0" euler="15 5 -10" childclass="link_geom">
        <site name="base_marker" pos="0.02 0.03 0.1" zaxis="0.3 -0.4 0.9"/>
        <body name="moving_link" pos="0.2 0 0.3">
          <frame pos="0.04 -0.06 0.05" axisangle="0 1 0 30">
            <joint name="pivot" type="hinge" pos="0.05 0.02 -0.01" axis="1 1 2" ref="23"/>
            <joint name="slider" type="slide" pos="0.02 -0.03 0.04" axis="0 2 -1" ref="0.12"/>
            <geom type="capsule" fromto="0 0 0 0.1 0.06 0.4" size="0.06" mass="0.9"/>
            <geom type="box" size="0.03 0.05 0.08" pos="0.02 -0.11 0.25" euler="0 25 15"/>
          </frame>
          <body name="payload" pos="0.11 0.07 0.42" euler="8 14 -21">
            <inertial pos="0.08 -0.02 0.05" mass="0.7" fullinertia="0.02 0.025 0.03 0.003 -0.002 0.001"/>
            <frame pos="0.02 0.01 0.06" axisangle="1 0 0 -18">
              <site name="tool" type="capsule" fromto="0 0 0 0.02 0.06 0.14" size="0.009"/>
            </frame>
          </body>
        </body>
      </body>
    </frame>
  </worldbody>
</mujoco>
"""


def pose(position, rotation=None):
    transform = np.eye(4)
    transform[:3, 3] = position
    if rotation is not None:
        transform[:3, :3] = rotation.reshape(3, 3)
    return transform


def rotation(quaternion):
    entries = np.empty(9)
    mujoco.mju_quat2Mat(entries, quaternion)
    return entries.reshape(3, 3)


def endpoint(model, label):
    if label.startswith("site:"):
        site = model.site(label[5:]).id
        return int(model.site_bodyid[site]), site
    return model.body(label).id, None


def world_pose(data, body, site):
    if site is None:
        return pose(data.xpos[body], data.xmat[body])
    return pose(data.site_xpos[site], data.site_xmat[site])


def numbers(values):
    return " ".join(format(float(value), ".17g") for value in values)


def add_origin(element, transform):
    r = transform[:3, :3]
    horizontal = math.hypot(r[0, 0], r[1, 0])
    pitch = math.atan2(-r[2, 0], horizontal)
    if horizontal < 1e-12:
        roll, yaw = 0, math.atan2(-r[0, 1], r[1, 1])
    else:
        roll, yaw = math.atan2(r[2, 1], r[2, 2]), math.atan2(r[1, 0], r[0, 0])
    ET.SubElement(element, "origin", xyz=numbers(transform[:3, 3]), rpy=numbers([roll, pitch, yaw]))


def add_joint(robot, name, kind, parent, child, transform):
    joint = ET.SubElement(robot, "joint", name=name, type=kind)
    ET.SubElement(joint, "parent", link=parent)
    ET.SubElement(joint, "child", link=child)
    add_origin(joint, transform)
    return joint


def add_body_link(robot, model, body):
    name = model.body(body).name or f"body_{body}"
    link = ET.SubElement(robot, "link", name=name)
    if model.body_mass[body] > 0:
        inertial = ET.SubElement(link, "inertial")
        add_origin(inertial, pose(model.body_ipos[body]))
        ET.SubElement(inertial, "mass", value=numbers([model.body_mass[body]]))
        axes = rotation(model.body_iquat[body])
        tensor = (axes * model.body_inertia[body]) @ axes.T
        keys = ("ixx", "ixy", "ixz", "iyy", "iyz", "izz")
        ET.SubElement(inertial, "inertia", dict(zip(keys, numbers(tensor[np.triu_indices(3)]).split())))
    return name


def export_urdf(model, data, label, base, flange, chain, output):
    base_name = model.site(base[1]).name if base[1] is not None else model.body(base[0]).name
    robot = ET.Element("robot", name=label)
    ET.SubElement(robot, "link", name=base_name)
    parent = base_name
    offset = np.linalg.solve(world_pose(data, *base), world_pose(data, base[0], None))
    joints = []
    for body in chain:
        offset = offset @ pose(model.body_pos[body], rotation(model.body_quat[body]))
        first = model.body_jntadr[body]
        for index in range(first, first + model.body_jntnum[body]):
            joint_type = model.jnt_type[index]
            if joint_type == mujoco.mjtJoint.mjJNT_HINGE:
                kind = "revolute" if model.jnt_limited[index] else "continuous"
            elif joint_type == mujoco.mjtJoint.mjJNT_SLIDE:
                kind = "prismatic"
            else:
                raise ValueError("Expected only hinge and slide joints in the test chain")
            offset_to_anchor = pose(model.jnt_pos[index])
            child = f"joint_{index}_frame"
            ET.SubElement(robot, "link", name=child)
            joint = add_joint(robot, model.joint(index).name, kind, parent, child, offset @ offset_to_anchor)
            ET.SubElement(joint, "axis", xyz=numbers(model.jnt_axis[index]))
            if kind != "continuous":
                reference = model.qpos0[model.jnt_qposadr[index]]
                low, high = model.jnt_range[index] - reference
                # URDF requires effort and velocity limits. The readers do not use them.
                ET.SubElement(joint, "limit", lower=numbers([low]), upper=numbers([high]), effort="0", velocity="0")
            joints.append(index)
            parent, offset = child, pose(-model.jnt_pos[index])
        child = add_body_link(robot, model, body)
        add_joint(robot, f"body_{body}_fixed", "fixed", parent, child, offset)
        parent, offset = child, np.eye(4)
    if flange[1] is not None:
        tip = model.site(flange[1]).name
        ET.SubElement(robot, "link", name=tip)
        site_pose = pose(model.site_pos[flange[1]], rotation(model.site_quat[flange[1]]))
        add_joint(robot, f"site_{flange[1]}_fixed", "fixed", parent, tip, offset @ site_pose)
        parent = tip
    ET.indent(robot, space="  ")
    ET.ElementTree(robot).write(output, encoding="utf-8", xml_declaration=True)
    return base_name, parent, np.array(joints, dtype=int)


def reference_samples(model, data, base, flange, joints):
    qpos = model.jnt_qposadr[joints]
    dofs = model.jnt_dofadr[joints]
    armature = np.diag(model.dof_armature[dofs])
    rng = np.random.default_rng(2026)
    scale = np.where(model.jnt_type[joints] == mujoco.mjtJoint.mjJNT_SLIDE, 0.12, 0.6)
    samples = []
    for i in range(13):
        inputs = np.zeros((3, len(joints))) if i == 0 else rng.uniform(-1, 1, (3, len(joints)))
        q, qd, qdd = inputs * scale
        data.qpos[:] = model.qpos0
        data.qpos[qpos] += q
        data.qvel[:] = 0
        data.qvel[dofs] = qd
        mujoco.mj_forward(model, data)
        relative_pose = np.linalg.solve(world_pose(data, *base), world_pose(data, *flange))
        mass = np.empty((len(dofs), len(dofs)))
        basis, column = np.zeros(model.nv), np.empty(model.nv)
        for j, dof in enumerate(dofs):
            basis[dof] = 1
            mujoco.mj_mulM(model, data, column, basis)
            mass[:, j] = column[dofs]
            basis[dof] = 0
        mass -= armature
        force = mass @ qdd + data.qfrc_bias[dofs]
        samples.append({
            "q": q.tolist(), "qd": qd.tolist(), "qdd": qdd.tolist(),
            "pose": relative_pose.ravel(order="F").tolist(),
            "mass": mass.ravel(order="F").tolist(), "force": force.tolist(),
        })
    return samples


def make_reference(label, xml, base_label, flange_label, output):
    model = mujoco.MjModel.from_xml_path(str(xml))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    base = endpoint(model, base_label)
    flange = endpoint(model, flange_label)
    chain = []
    body = flange[0]
    while body != base[0]:
        if body == 0:
            raise ValueError("Expected the flange to descend from the base")
        chain.append(body)
        body = model.body_parentid[body]
    urdf = output / f"{label}.urdf"
    base_name, tip, joints = export_urdf(model, data, label, base, flange, reversed(chain), urdf)
    gravity = world_pose(data, *base)[:3, :3].T @ model.opt.gravity
    return {
        "label": label, "mjcf": str(xml), "urdf": str(urdf),
        "base": base_label, "flange": flange_label, "urdf_base": base_name, "urdf_flange": tip,
        "joint_names": [model.joint(i).name for i in joints],
        "anchors": model.jnt_pos[joints].ravel().tolist(), "gravity": gravity.tolist(),
        "samples": reference_samples(model, data, base, flange, joints),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    output = parser.parse_args().output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    repo = Path(__file__).resolve().parents[1]
    fixture = output / "joint_frames.xml"
    fixture.write_text(FRAME_FIXTURE, encoding="utf-8")
    cases = (
        ("ur5e", repo / "mj_models/universal_robots_ur5e/scene.xml", "world", "site:attachment_site"),
        ("iiwa14", repo / "mj_models/kuka_iiwa_14/scene.xml", "world", "site:attachment_site"),
        ("cart_pole", repo / "simulink/mj_IP.xml", "base_link", "rod"),
        ("joint_frames", fixture, "site:base_marker", "site:tool"),
    )
    reference = [make_reference(*case, output) for case in cases]
    (output / "reference.json").write_text(json.dumps(reference, allow_nan=False), encoding="utf-8")
    print(f"MuJoCo {mujoco.__version__}: exported {len(reference)} URDF chains.")


if __name__ == "__main__":
    main()
