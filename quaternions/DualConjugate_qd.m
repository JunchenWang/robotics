function qd_bar = DualConjugate_qd(qd)
qd_bar = [Conjugate_q(qd(1 : 4))
          -Conjugate_q(qd(5 : 8))];