def f(x, y):
    return x + y


def rk4_step(x, y, h):
    k1 = h * f(x, y)
    k2 = h * f(x + h/2, y + k1/2)
    k3 = h * f(x + h/2, y + k2/2)
    k4 = h * f(x + h, y + k3)

    return y + (k1 + 2*k2 + 2*k3 + k4) / 6


x0 = 0.0
y0 = 1.0
h = 0.1
n = 5

x = [x0]
y = [y0]

# Calculate the first three steps using RK4
for i in range(3):
    y_next = rk4_step(x[i], y[i], h)
    x.append(x[i] + h)
    y.append(y_next)

# Apply the four-step Adams-Bashforth formula
for i in range(3, n):
    y_next = y[i] + (h / 24) * (
        55 * f(x[i], y[i])
        - 59 * f(x[i-1], y[i-1])
        + 37 * f(x[i-2], y[i-2])
        - 9 * f(x[i-3], y[i-3])
    )

    x.append(x[i] + h)
    y.append(y_next)


print("X\tY")

for i in range(n + 1):
    print(f"{x[i]:.1f}\t{y[i]:.5f}")