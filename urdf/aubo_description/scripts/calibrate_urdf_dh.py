#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
用途
----
- 从 URDF 中提取各关节在 q=0 时的 Modified DH 参数。
- 通过 pyaubo_sdk RPC 获取运动学校准补偿值 (da, dalpha, dd, dtheta)。
- 将补偿后的结果回写为新的 URDF 文件，默认生成“校准版”文件而不是覆盖原文件。

示例
----
python3 aubo_description/scripts/calibrate_urdf_dh.py --robot-model aubo_i5 --robot-ip 172.16.31.2 --output-path aubo_description/urdf/aubo_i5_calibrated.urdf

python3 scripts/calibrate_urdf_dh.py \
  --robot-model aubo_i10H \
  --robot-ip 192.168.15.128 \
  --temperature 20

python3 scripts/calibrate_urdf_dh.py \
  --urdf-in aubo_i10H \
  --output-path urdf/aubo_i10H_calibrated.urdf \
  --robot-ip 192.168.15.128
"""

from __future__ import annotations

import argparse
import ast
import json
import math
import os
import re
import sys
from collections import defaultdict
from typing import Dict, List, Tuple
import xml.etree.ElementTree as ET

try:
    import numpy as np
    NUMPY_IMPORT_ERROR = None
except Exception as exc:
    np = None
    NUMPY_IMPORT_ERROR = exc

try:
    import pyaubo_sdk  # optional at import-time
    PYAUBO_SDK_IMPORT_ERROR = None
except Exception:
    pyaubo_sdk = None  # allow --help and offline checks without the SDK
    PYAUBO_SDK_IMPORT_ERROR = sys.exc_info()[1]


def format_dependency_error(module_name: str, import_error: Exception | None, hint: str) -> str:
    parts = ["缺少运行依赖 `{}`。".format(module_name)]
    if import_error is not None:
        parts.append("导入错误: {}".format(import_error))
    parts.append(hint)
    return " ".join(parts)


def check_runtime_dependencies():
    errors = []
    if np is None:
        errors.append(
            format_dependency_error(
                "numpy",
                NUMPY_IMPORT_ERROR,
                "当前环境无法导入 numpy，请先确认运行环境依赖是否齐全。",
            )
        )
    if pyaubo_sdk is None:
        errors.append(
            format_dependency_error(
                "pyaubo_sdk",
                PYAUBO_SDK_IMPORT_ERROR,
                "当前环境无法导入 pyaubo_sdk，请先确认 AUBO Python SDK 是否已正确配置到当前终端环境。",
            )
        )
    if errors:
        raise RuntimeError("依赖检查失败:\n- " + "\n- ".join(errors))


def rpy_to_matrix(roll: float, pitch: float, yaw: float) -> np.ndarray:
    """URDF fixed-axis order: R = Rz(yaw) * Ry(pitch) * Rx(roll)."""
    cr, sr = math.cos(roll), math.sin(roll)
    cp, sp = math.cos(pitch), math.sin(pitch)
    cy, sy = math.cos(yaw), math.sin(yaw)
    rz = np.array([[cy, -sy, 0], [sy, cy, 0], [0, 0, 1]], dtype=float)
    ry = np.array([[cp, 0, sp], [0, 1, 0], [-sp, 0, cp]], dtype=float)
    rx = np.array([[1, 0, 0], [0, cr, -sr], [0, sr, cr]], dtype=float)
    return rz.dot(ry).dot(rx)


def matrix_to_rpy(rotation: np.ndarray) -> Tuple[float, float, float]:
    """Convert a rotation matrix to URDF fixed-axis RPY."""
    sy = -rotation[2, 0]
    cy = math.sqrt(max(0.0, 1.0 - sy * sy))
    singular = cy < 1e-8
    if not singular:
        pitch = math.asin(sy)
        roll = math.atan2(rotation[2, 1], rotation[2, 2])
        yaw = math.atan2(rotation[1, 0], rotation[0, 0])
    else:
        pitch = math.asin(sy)
        roll = 0.0
        yaw = math.atan2(-rotation[0, 1], rotation[1, 1])
    return roll, pitch, yaw


def rot_x(alpha: float) -> np.ndarray:
    ca, sa = math.cos(alpha), math.sin(alpha)
    return np.array([[1, 0, 0], [0, ca, -sa], [0, sa, ca]], dtype=float)


def rot_z(theta: float) -> np.ndarray:
    ct, st = math.cos(theta), math.sin(theta)
    return np.array([[ct, -st, 0], [st, ct, 0], [0, 0, 1]], dtype=float)


def parse_urdf(path: str):
    tree = ET.parse(path)
    root = tree.getroot()
    return tree, root


def get_joints(root: ET.Element) -> List[Dict]:
    joints = []
    for joint in root.findall("joint"):
        joint_type = joint.get("type", "")
        if joint_type not in ("revolute", "prismatic", "continuous", "fixed"):
            continue
        parent_el = joint.find("parent")
        child_el = joint.find("child")
        if parent_el is None or child_el is None:
            continue
        origin = joint.find("origin")
        if origin is None:
            xyz = [0.0, 0.0, 0.0]
            rpy = [0.0, 0.0, 0.0]
        else:
            xyz = [float(value) for value in origin.get("xyz", "0 0 0").split()]
            rpy = [float(value) for value in origin.get("rpy", "0 0 0").split()]
        joints.append(
            {
                "elem": joint,
                "name": joint.get("name"),
                "type": joint_type,
                "parent": parent_el.get("link"),
                "child": child_el.get("link"),
                "xyz": np.array(xyz, dtype=float),
                "rpy": np.array(rpy, dtype=float),
            }
        )
    return joints


def order_chain(joints: List[Dict], base_link: str | None = None):
    parent_to_joint = defaultdict(list)
    link_parents = {}
    links = set()
    for joint in joints:
        parent_to_joint[joint["parent"]].append(joint)
        link_parents[joint["child"]] = joint["parent"]
        links.add(joint["parent"])
        links.add(joint["child"])

    if base_link is None:
        candidates = [link for link in links if link not in link_parents]
        if not candidates:
            raise ValueError("无法自动识别 base_link，请显式传入 --base-link。")
        base_link = candidates[0]

    ordered = []
    current = base_link
    seen = {current}
    while True:
        candidates = parent_to_joint.get(current, [])
        if not candidates:
            break
        ordered_joint = sorted(candidates, key=lambda item: (item["type"] == "fixed",))[0]
        ordered.append(ordered_joint)
        current = ordered_joint["child"]
        if current in seen:
            break
        seen.add(current)
    return ordered, base_link


def modified_decompose(rotation: np.ndarray, position: np.ndarray):
    alpha = math.atan2(rotation[2, 1], rotation[2, 2])
    rz_component = rot_x(-alpha).dot(rotation)
    theta0 = math.atan2(rz_component[1, 0], rz_component[0, 0])
    a_val = float(position[0])
    py, pz = float(position[1]), float(position[2])
    d_val = math.hypot(py, pz)
    return a_val, alpha, d_val, theta0


def modified_compose(a_val: float, alpha: float, d_val: float, theta0: float):
    rotation = rot_x(alpha).dot(rot_z(theta0))
    position = np.array(
        [a_val, -d_val * math.sin(alpha), d_val * math.cos(alpha)], dtype=float
    )
    return rotation, position


def extract_transforms(joints_ordered: List[Dict]):
    transforms = []
    for joint in joints_ordered:
        if joint["type"] == "fixed":
            continue
        transforms.append(
            (rpy_to_matrix(*joint["rpy"]), joint["xyz"].astype(float), joint["name"])
        )
    return transforms


def load_deltas_str(delta_input, joints_ordered: List[Dict]) -> Dict[str, Dict[str, float]]:
    """Parse the SDK compensation payload."""
    if isinstance(delta_input, dict):
        data = delta_input
    elif isinstance(delta_input, str):
        text = delta_input.strip()
        if not text:
            raise ValueError("校准补偿为空。")
        try:
            data = json.loads(text)
        except Exception:
            data = ast.literal_eval(text)
    else:
        raise TypeError("delta_input 必须是 dict 或 str。")

    if not isinstance(data, dict):
        raise ValueError(
            "校准补偿格式错误，应类似 {'a': [...], 'alpha': [...], 'd': [...], 'theta': [...]}。"
        )

    joint_names = [joint["name"] for joint in joints_ordered if joint["type"] != "fixed"]
    param_map = {"a": "da", "alpha": "dalpha", "d": "dd", "theta": "dtheta"}

    def to_float_list(value):
        if value is None:
            return []
        sequence = value if isinstance(value, (list, tuple)) else [value]
        return [float(item) for item in sequence]

    arrays = {key: to_float_list(value) for key, value in data.items() if key in param_map}
    if not arrays:
        raise ValueError("校准补偿中缺少 a/alpha/d/theta。")

    deltas = {}
    for index, joint_name in enumerate(joint_names):
        delta = {}
        for source_key, target_key in param_map.items():
            values = arrays.get(source_key)
            if values is not None and index < len(values):
                delta[target_key] = values[index]
        deltas[joint_name] = delta
    return deltas


def load_deltas_txt_matrix(path: str, joints_ordered: List[Dict]) -> Dict[str, Dict[str, float]]:
    """Optional MATLAB-like text parser kept for local debugging."""
    joint_names = [joint["name"] for joint in joints_ordered if joint["type"] != "fixed"]
    matlab_arrays = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith(";"):
                continue
            match = re.match(r"delta_(\w+)\s*=\s*\[([^\]]+)\];", line)
            if match:
                param_name = match.group(1)
                values = [float(item.strip()) for item in match.group(2).split(",")]
                matlab_arrays[param_name] = values

    if not matlab_arrays:
        raise ValueError("未在文本文件中读取到校准补偿数组。")

    param_map = {"a": "da", "alpha": "dalpha", "d": "dd", "theta": "dtheta"}
    deltas = {}
    for index, joint_name in enumerate(joint_names):
        delta = {}
        for matlab_param, values in matlab_arrays.items():
            if matlab_param in param_map and index < len(values):
                delta[param_map[matlab_param]] = values[index]
        deltas[joint_name] = delta
    return deltas


def print_modified_dh_params(transforms):
    print("URDF 当前 Modified DH 参数 (q=0):")
    for rotation, position, name in transforms:
        a_val, alpha, d_val, theta0 = modified_decompose(rotation, position)
        print(
            "{}: a={:.9f}, alpha={:.9f}, d={:.9f}, theta0={:.9f}".format(
                name, a_val, alpha, d_val, theta0
            )
        )


def print_modified_dh_after_deltas(transforms, deltas: Dict[str, Dict[str, float]]):
    print("\n应用校准补偿后的 Modified DH 参数 (m/rad):")
    for rotation, position, name in transforms:
        a_val, alpha, d_val, theta0 = modified_decompose(rotation, position)
        delta = deltas.get(name, {})
        a_new = a_val + float(delta.get("da", 0.0))
        alpha_new = alpha + float(delta.get("dalpha", 0.0))
        d_new = d_val + float(delta.get("dd", 0.0))
        theta_new = theta0 + float(delta.get("dtheta", 0.0))
        print(
            "{}: a={:.9f}, alpha={:.9f}, d={:.9f}, theta0={:.9f}".format(
                name, a_new, alpha_new, d_new, theta_new
            )
        )

    print("\n应用校准补偿后的 Modified DH 参数 (mm/deg):")
    for rotation, position, name in transforms:
        a_val, alpha, d_val, theta0 = modified_decompose(rotation, position)
        delta = deltas.get(name, {})
        a_new = a_val + float(delta.get("da", 0.0))
        alpha_new = alpha + float(delta.get("dalpha", 0.0))
        d_new = d_val + float(delta.get("dd", 0.0))
        theta_new = theta0 + float(delta.get("dtheta", 0.0))
        print(
            "{}: a={:.3f}, alpha={:.3f}, d={:.3f}, theta0={:.3f}".format(
                name,
                a_new * 1000.0,
                math.degrees(alpha_new),
                d_new * 1000.0,
                math.degrees(theta_new),
            )
        )


def ensure_parent_dir(path: str):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def infer_package_root() -> str:
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def resolve_urdf_candidate(path_or_name: str, package_root: str) -> str:
    value = path_or_name.strip()
    if not value:
        return ""

    if os.path.isabs(value) or os.path.sep in value:
        return os.path.abspath(value)

    if value.endswith(".urdf"):
        return os.path.join(package_root, "urdf", value)

    return os.path.join(package_root, "urdf", "{}.urdf".format(value))


def resolve_input_urdf(args, package_root: str) -> str:
    if args.urdf_in:
        return resolve_urdf_candidate(args.urdf_in, package_root)
    if args.robot_model:
        return resolve_urdf_candidate(args.robot_model, package_root)
    raise ValueError("请至少提供 --robot-model 或 --urdf-in，可直接传入机器人名称。")


def build_default_output_path(input_urdf: str, package_root: str, suffix: str) -> str:
    input_abs = os.path.abspath(input_urdf)
    input_name = os.path.splitext(os.path.basename(input_abs))[0]
    out_dir = os.path.join(package_root, "urdf")
    return os.path.join(out_dir, "{}{}.urdf".format(input_name, suffix))


def indent_xml(tree: ET.ElementTree):
    if hasattr(ET, "indent"):
        ET.indent(tree, space="  ")


def fetch_compensation_from_rpc(args):
    rpc = pyaubo_sdk.RpcClient()
    rpc.connect(args.robot_ip, args.robot_port)
    if not rpc.hasConnected():
        raise RuntimeError(
            "RPC 连接失败: {}:{}".format(args.robot_ip, args.robot_port)
        )
    print("已连接到机器人 RPC: {}:{}".format(args.robot_ip, args.robot_port))

    rpc.login(args.user, args.password)
    if not rpc.hasLogined():
        raise RuntimeError("RPC 登录失败: user={}".format(args.user))
    print("RPC 登录成功。")

    robot_names = rpc.getRobotNames()
    if not robot_names:
        raise RuntimeError("RPC 未返回机器人名称。")
    robot_if = rpc.getRobotInterface(robot_names[0])
    cfg = robot_if.getRobotConfig()

    robot_type = cfg.getRobotType()
    robot_subtype = cfg.getRobotSubType()
    print("机器人 type:", robot_type)
    print("机器人 subtype:", robot_subtype)

    compensation = cfg.getKinematicsCompensate(args.temperature)
    print("原始校准补偿数据:", compensation)
    return compensation


def apply_calibration_to_tree(tree, root, joints_ordered, deltas):
    transforms = extract_transforms(joints_ordered)
    name_to_joint = {joint["name"]: joint for joint in joints_ordered}

    for rotation, position, name in transforms:
        a_val, alpha, d_val, theta0 = modified_decompose(rotation, position)
        delta = deltas.get(name, {})
        a_new = a_val + float(delta.get("da", 0.0))
        alpha_new = alpha + float(delta.get("dalpha", 0.0))
        d_new = d_val + float(delta.get("dd", 0.0))
        theta_new = theta0 + float(delta.get("dtheta", 0.0))

        rotation_new, position_new = modified_compose(a_new, alpha_new, d_new, theta_new)
        roll, pitch, yaw = matrix_to_rpy(rotation_new)

        joint = name_to_joint.get(name)
        if joint is None or joint["type"] == "fixed":
            continue
        origin_el = joint["elem"].find("origin")
        if origin_el is None:
            origin_el = ET.SubElement(joint["elem"], "origin")
        origin_el.set(
            "xyz",
            "{:.9f} {:.9f} {:.9f}".format(
                position_new[0], position_new[1], position_new[2]
            ),
        )
        origin_el.set("rpy", "{:.16f} {:.16f} {:.16f}".format(roll, pitch, yaw))

    root.set("name", root.get("name", "robot"))
    indent_xml(tree)


def parse_args():
    parser = argparse.ArgumentParser(
        description="读取机器人校准补偿并生成新的校准版 URDF 文件。"
    )
    parser.add_argument("--robot-model", type=str, default="", help="机器人型号，例如 aubo_i10H。")
    parser.add_argument(
        "--urdf-in",
        type=str,
        default="",
        help="输入 URDF 路径或机器人名称，例如 urdf/aubo_i10H.urdf 或 aubo_i10H。",
    )
    parser.add_argument(
        "--output-path",
        type=str,
        default="",
        help="输出 URDF 路径。默认生成到 aubo_description/urdf/<robot>_calibrated.urdf。",
    )
    parser.add_argument(
        "--output-suffix",
        type=str,
        default="_calibrated",
        help="默认输出文件后缀，默认值为 _calibrated。",
    )
    parser.add_argument("--base-link", type=str, default="", help="可选，显式指定 base link 名称。")
    parser.add_argument("--temperature", type=float, default=20.0, help="获取校准补偿时使用的温度。")
    parser.add_argument("--robot-ip", type=str, required=True, help="机器人 RPC IP，必须显式传入。")
    parser.add_argument("--robot-port", type=int, default=30004, help="机器人 RPC 端口。")
    parser.add_argument("--user", type=str, default="aubo", help="机器人 RPC 用户名。")
    parser.add_argument("--password", type=str, default="123456", help="机器人 RPC 密码。")
    parser.add_argument(
        "--skip-dependency-check",
        action="store_true",
        help="跳过启动时的依赖检查，不推荐，通常仅用于调试。",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="允许覆盖已存在的输出文件；默认禁止覆盖，避免误改原始文件。",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if not args.skip_dependency_check:
        check_runtime_dependencies()

    package_root = infer_package_root()
    input_urdf = resolve_input_urdf(args, package_root)
    if not os.path.isfile(input_urdf):
        raise FileNotFoundError("输入 URDF 不存在: {}".format(input_urdf))

    output_path = args.output_path.strip()
    if output_path:
        output_path = os.path.abspath(output_path)
    else:
        output_path = build_default_output_path(input_urdf, package_root, args.output_suffix)

    if os.path.abspath(output_path) == os.path.abspath(input_urdf):
        raise ValueError("输出文件不能与输入 URDF 相同，请生成新的校准文件。")
    if os.path.exists(output_path) and not args.force:
        raise FileExistsError(
            "输出文件已存在: {}。如需覆盖，请显式传入 --force。".format(output_path)
        )

    tree, root = parse_urdf(input_urdf)
    print("已加载 URDF:", input_urdf)

    joints = get_joints(root)
    joints_ordered, base_link = order_chain(joints, args.base_link or None)
    print("识别到 base_link:", base_link)

    transforms = extract_transforms(joints_ordered)
    print_modified_dh_params(transforms)

    compensation = fetch_compensation_from_rpc(args)
    deltas = load_deltas_str(compensation, joints_ordered) if compensation else {}
    print_modified_dh_after_deltas(transforms, deltas)

    apply_calibration_to_tree(tree, root, joints_ordered, deltas)
    ensure_parent_dir(output_path)
    tree.write(output_path, encoding="utf-8", xml_declaration=True)
    print("\n已生成校准版 URDF:", output_path)
    print("请重新执行 `colcon build --packages-select aubo_description`，使新生成的 URDF 被安装到工作空间。")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("Error:", exc)
        sys.exit(1)
