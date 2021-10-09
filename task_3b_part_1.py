with open('Task_3b.txt') as inf:
    nx = int(inf.readline().rstrip().replace('nx = ', ''))  # Число точек, в которых производились измерения
    ny = int(inf.readline().rstrip().replace('ny = ', ''))  # Количество измерений в каждой точке
    x = [float(i) for i in inf.readline().rstrip().replace('x= ', '').split(' ')]  # Сами точки
    y = dict()
    for xi in x:
        line = inf.readline().rstrip().split('=')[1]
        y[xi] = [float(i) for i in line.split(' ')]  # Все измерения (n штук) в соответствующей точке

with open('x.txt', 'w') as ouf:
    for xi in x:
        ouf.write(str(xi) + ' ')

with open('y.txt', 'w') as ouf:
    for xi in x:
        for yi in y[xi]:
            ouf.write(str(yi) + ' ')
        ouf.write('\n')
