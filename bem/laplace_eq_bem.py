from typing import List

import numpy as np
from numpy.typing import NDArray

from bem.func_settings import make_points, normal_vector
from bem.utils import EdgePoints, exact_u, integral, solve

np.seterr(divide="raise")


class Problem(EdgePoints):
    """問題設定クラス"""

    def __init__(self, div_num) -> None:
        self.div_num = div_num
        # 点列を生成. これでU, Wの要素も計算できる
        super().__init__(make_points(div_num))
        # 境界条件を生成
        self.bc = self.make_boundary_condition()
        # 影響係数行列を計算
        self.U_mat = np.array(
            [
                [self.calc_U(m, n) for n in range(self.div_num)]
                for m in range(self.div_num)
            ]
        )
        self.W_mat = np.array(
            [
                [self.calc_W(m, n) for n in range(self.div_num)]
                for m in range(self.div_num)
            ]
        )
        # 共役境界量を計算
        self.conj_bd_values = self.calc_conjugate_boundary_values()

    def make_boundary_condition(self) -> NDArray[np.float64]:
        """
        点列を受け取りそれらに対する境界条件を作成.
        pointsのshapeは(div_num, 2).
        返り値のshapeは(div_num,).
        """
        # とりあえずDirichlet条件
        return np.array([exact_u(x) for x in self.points])

    def calc_conjugate_boundary_values(self) -> NDArray[np.float64]:
        """
        影響係数行列と既知の境界条件から未知の共役な境界量を求める.
        返り値のshapeは(div_num,).
        """
        solution = solve(self.U_mat, np.dot(self.W_mat, self.bc))
        return solution

    def interior_point_calc(self, x: NDArray[np.float64]):
        """内点xの値の解を計算する"""
        y = self.points
        # 内点計算に使う基本解の値の列を計算
        fund_val = -np.log(np.linalg.norm(x - y, axis=1)) / 2.0 / np.pi
        # 基本解の法線微分の列も計算
        n_vecs = np.array([normal_vector(point) for point in y])
        fund_dv_val = (
            np.einsum("ij,ij->i", x - y, n_vecs)
            / 2.0
            / np.pi
            / np.linalg.norm(x - y, axis=1)
        )
        # 積分するための配列
        # とりあえず半径1の円周（あとでself.pointsから積分を構成できるように変更する）
        circ_x = 2 * np.pi * np.arange(self.div_num + 1) / self.div_num
        result = 0
        # 一重層ポテンシャル; 周期的な配列を作成
        circ_y = np.append(
            fund_val * self.conj_bd_values, fund_val[0] * self.conj_bd_values[0]
        )
        result += integral(circ_y, circ_x)
        # 二重層ポテンシャル
        circ_y = np.append(fund_dv_val * self.bc, fund_dv_val[0] * self.bc[0])
        result -= integral(circ_y, circ_x)

        return result


def main(div_num):
    """境界要素法を使用してラプラス問題を解く"""
    # 影響係数行列を作成. ここに点列情報も入っている.
    problem = Problem(div_num)
    # 計算する内点リストを作成.
    r_div = 8
    calc_points = np.array(
        [
            [
                r_num / r_div * np.cos(2 * np.pi * th_num / div_num),
                r_num / r_div * np.sin(2 * np.pi * th_num / div_num),
            ]
            for r_num in range(1, r_div)
            for th_num in range(div_num)
        ]
    )
    # 計算する内点に原点を追加
    calc_points = np.append(np.array([[0.0, 0.0]]), calc_points, axis=0)
    # 各点について値を計算
    result: List[float] = []
    for x in calc_points:
        result.append(problem.interior_point_calc(x))
    result_ndarray = np.array(result)

    # 境界上の点についての情報を追加して返す
    point_list = np.append(calc_points, problem.points, axis=0)
    value_list = np.append(result_ndarray, problem.bc)

    return point_list, value_list
