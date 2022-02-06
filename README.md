# 境界要素法（BEM）でラプラス方程式の境界値問題を解く

## 問題設定
<img src="https://latex.codecogs.com/svg.image?\Delta&space;u=0" title="\Delta u=0" />
という問題を
<img src="https://latex.codecogs.com/svg.image?\bar{u}=x^3-3xy^2" title="\bar{u}=x^3-3xy^2" />
というDirichlet境界条件のもとで解く.
境界は半径1, 原点中心の単位円.

## 実行
1. `pipenv install`
1. `pipenv shell`
1. `python main.py`

により, 解いた結果がプロットされる.

## その他
領域や関数を変更したい場合は`bem/func_settings.py`の中を変更する
~~（いまのところ`bem/laplace_eq_bem.py`の`interior_point_calc`の`circ_x`を定義する部分も変更しないといけない）~~.
後述の領域刻み数が大きい場合は線積分は概ねうまく自動で実行できるが,
正確さが求められる場合は`bem/func_settings.py`の`linear_integral_x`で積分領域の数値列を定義するなどし, `bem/laplace_eq_bem.py`の`Problem`クラスのコンストラクタで`self.circ_x`に格納しておく.

境界条件を変更する場合は, `bem/laplace_eq_bem.py`の`make_boundary_condition`を変更.
計算する内点の列は`bem/laplace_eq_bem.py`の`main`で定義しているので, ここを変更する.
この`main`は領域刻み数を引数として受け取り, 内点の配列と各点での値の配列をタプルで返す.