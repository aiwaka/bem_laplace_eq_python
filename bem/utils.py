from dataclasses import dataclass

import numpy as np
from numpy.linalg import LinAlgError
from numpy.typing import NDArray
from scipy import integrate

from bem.func_settings import grad_u, normal_vector, u


@dataclass
class ComponentValues:
    x_m: NDArray
    x1: NDArray
    x2: NDArray
    lX1: float
    lX2: float
    lY1: float
    lY2: float
    r1: float
    r2: float
    h: float
    theta: float


class EdgePoints:
    def __init__(self, points: NDArray[np.float64]) -> None:
        """pointsは二次元空間点の配列, p_numは点数"""
        self.points = points
        self.p_num = points.shape[0]

    def calc_component_values(self, m, n) -> ComponentValues:
        """pointsと添字m,nから影響係数に使う値を計算"""
        x_m = (self.points[m % self.p_num] + self.points[(m + 1) % self.p_num]) / 2
        # n番目の区間の端点
        x1 = self.points[n % self.p_num]
        x2 = self.points[(n + 1) % self.p_num]

        # 概ね小林本の表式に合わせ, XやYを導入している.
        h = np.linalg.norm(x1 - x2, axis=0)  # 区間の長さ

        t_vec = (x2 - x1) / h
        n_vec = np.array([x1[1] - x2[1], x2[0] - x1[0]]) / h

        lX1 = np.dot(x_m - x1, t_vec)
        lX2 = np.dot(x_m - x2, t_vec)
        r1 = np.linalg.norm(x_m - x1, axis=0)  # 型アノテーションをAnyにするためaxisを指定している
        r2 = np.linalg.norm(x_m - x2, axis=0)
        lY1 = np.dot(x_m - x1, n_vec)
        lY2 = np.dot(x_m - x2, n_vec)
        theta = np.arctan2(lY2, lX2) - np.arctan2(lY1, lX1)

        return ComponentValues(x_m, x1, x2, lX1, lX2, lY1, lY2, r1, r2, h, theta)

    def calc_U(self, m, n):
        """添字m, nに対する影響係数Uを計算"""
        val = self.calc_component_values(m, n)
        if m == n:
            # m=nの場合は例外的な計算
            return (1.0 - np.log(val.h / 2)) * val.h / 2 / np.pi

        return (
            (
                val.lX2 * np.log(val.r2)
                - val.lX1 * np.log(val.r1)
                + val.h
                - val.lY1 * val.theta
            )
            / 2.0
            / np.pi
        )

    def calc_W(self, m, n):
        """m, nに対するWを計算"""
        if m == n:
            return 0.5
        val = self.calc_component_values(m, n)

        return val.theta / 2.0 / np.pi


def exact_u(x):
    """
    二次元座標xを受け取り, uの値を返す
    """
    return u(x)


def exact_u_normal_dv(x):
    """二次元座標xを受け取り, 法線微分の値を返す"""
    return np.dot(grad_u(x), normal_vector(x))


def fund_gamma(x, y):
    """基本解の値を返す"""
    return -np.log(np.linalg.norm(x - y)) / 2.0 / np.pi


def fund_gamma_normal_dv(x, y):
    return np.dot(x - y, normal_vector(y)) / 2.0 / np.pi / np.linalg.norm(x - y)


def solve(A, b):
    """Ax=bの方程式を解く. Aがsingularの場合例外を出す."""
    try:
        x = np.linalg.solve(A, b)
    except LinAlgError as e:
        print(e)
        exit(1)
    return x


def integral(y, x):
    """xにより離散化されたyの値を使って数値積分を行う"""
    # 台形則の場合累積和の配列になるので最後の要素を結果として返す
    return integrate.cumulative_trapezoid(y, x)[-1]
