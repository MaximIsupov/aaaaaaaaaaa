def f1(x1, x2):
    return 2 * x1 - x2 + 1

def f2(x1, x2):
    return 8 * x2 - x1

def f3(x01, x02, x1, x2, alf):
    return -2 * x01 * (x1 - x01 * alf) - 8 * x02 * (x2 - x02 * alf)

x, e1, e2, n = {x1 = float(input()), x2 = float(input())}, float(input()), float(input()), int(input())
k = 0
ans = {}
prev = {}
while k < n:
    grad.x1 = f1(x[x1], x[x2])
    grad.x2 = f2(x[x1], x[x2])
    prev.x1 = x.x1
    prev.x2 = x.x2
    if max(x[x1], x[x2]) < e1:
        ans = x
        break
    else:
        a = 0
	eps = 0.001
	b = 1
	si = 0.5;
	iter = 0;
	while (a + b) / 2 > eps and iter < 1000:
		mid = (a + b) / 2;
        	if (f3(grad.x1, grad.x2, x.x1, x.x2, si - mid) > funct(grad.x1, grad.x2, x.x1, x.x2, si + mid)) {
			a = mid - si;
		}
		else {
			b = mid + si;
		}
		iter++;
	}
	minim = (a + b) / 2
	x.x1 -= grad.x1 * minim
	x.x2 -= grad.x2 * minim
	if k!= 0:
            max(
