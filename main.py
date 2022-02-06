import matplotlib.pyplot as plt

from bem import laplace_eq_bem


def main():
    points, z = laplace_eq_bem.main(128)
    x = points[:, 0]
    y = points[:, 1]
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.scatter(x, y, z)
    plt.show()


if __name__ == "__main__":
    main()
