function v2 = TransformByQ(q, v1)
v2 = ConcatenateQuaternions(q, q_vector(v1), Conjugate_q(q));
v2 = v2(1:3);