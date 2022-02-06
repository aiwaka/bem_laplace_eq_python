"""
用いる関数についての設定などを行う
"""

import numpy as np
from numpy.typing import NDArray


def u(x):
    """ここで考える関数uを定める."""
    return x[0] ** 3 - 3 * x[0] * (x[1] ** 2)


def grad_u(x):
    """uの勾配（x1,x2方向偏微分の組）を返す."""
    return np.array([3 * (x[0] ** 2 - x[1] ** 2), -6 * x[0] * x[1]])


def normal_vector(x):
    """境界領域上のxにおける法線ベクトルを返す"""
    theta = np.arctan2(x[1], x[0])
    return np.array([np.cos(theta), np.sin(theta)])


def make_points(div_num) -> NDArray[np.float64]:
    """
    ここで離散化された境界の点列を作る. 境界を変更したいならこれを変更する.
    返り値のshapeは(div_num, 2).
    """
    phases = 2 * np.pi * np.arange(div_num) / div_num
    return np.array([np.cos(phases), np.sin(phases)]).T
