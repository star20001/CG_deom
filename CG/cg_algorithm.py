import math

def draw_line(p_list,algorithm):
    '''
    :param p_list: list of list of int:[[x0,y0],[x1,y1]]
    :param algorithm:string,the way to draw
    :return:(list of list of int:[[x_0,y_0],[x_1,y_1]...]
    '''
    x0, y0 = p_list[0]
    x1, y1 = p_list[1]
    result = []
    if algorithm == 'Naive':
        if x0 == x1:
            for y in range(y0, y1 + 1):
                result.append((x0, y))
        else:
            if x0 > x1:
                x0, y0, x1, y1 = x1, y1, x0, y0
            k = (y1 - y0) / (x1 - x0)
            for x in range(x0, x1 + 1):
                result.append((x, int(y0 + k * (x - x0))))
    elif algorithm == 'DDA':
        length = max(abs(x1 - x0), abs(y1 - y0))
        if length == 0:
            delta_x = delta_y = 0
        else:
            delta_x = (x1 - x0) / length
            delta_y = (y1 - y0) / length

        x = x0 + 0.5
        y = y0 + 0.5

        i = 1
        while(i <= length):
            result.append((int(x), int(y)))
            x = x + delta_x
            y = y + delta_y
            i = i + 1

    elif algorithm == 'Bresenham':
        flag = False
        delta_x = abs(x1 - x0)
        delta_y = abs(y1 - y0)
        result.append((x0, y0))
        if delta_x == 0 and delta_y == 0:  # only a dot
            return result
        if(delta_x < delta_y):
            flag = True
            x0, y0 = y0, x0
            x1, y1 = y1, x1
            delta_x, delta_y = delta_y, delta_x

        if x1 - x0 > 0:
            tx = 1
        else:
            tx = -1
        if y1 - y0 > 0:
            ty = 1
        else:
            ty = -1

        x = x0
        y = y0
        e = 2 * delta_y - delta_x
        while x != x1:
            if e < 0:  # right
                e = e + 2 * delta_y
            else:  # top_right
                e = e + 2 * delta_y - 2 * delta_x
                y = y + ty
            x = x + tx
            if flag:
                result.append((y, x))
            else:
                result.append((x, y))

    return result


def draw_polygon(p_list,algorithm):
    '''
    :param p_list: (list of list of int:[[x0,y0],[x1,y1]...]
    :param algorithm:(string)
    :return:(list of list of int:[[x0,y0]...]
    '''
    result=[]
    for i in range(len(p_list)):
        line=draw_line([p_list[i-1],p_list[i]],algorithm)
        result+=line
    return result

def draw_polygon_gui(p_list,algorithm):
    '''
    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...])
    :param algorithm: (string) such as 'DDA' or 'Bresenham'
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...])
    '''
    result=[]
    for i in range(len(p_list)):
        line=draw_line([p_list[i-1],p_list[i]],algorithm)
        result+=line
    return result

def draw_ellipse(p_list):
    '''
    :param p_list: (list of list of int: [[x0, y0], [x1, y1]])
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...])
    '''
    x0, y0 = p_list[0]
    x1, y1 = p_list[1]
    xc = (x0 + x1) // 2
    yc = (y0 + y1) // 2  #center
    rx = abs((x1 - x0)) // 2
    ry = abs((y1 - y0)) // 2
    f = False
    if rx < ry:
        rx, ry = ry, rx
        f = True
    x, y = 0, ry
    ry_sq, rx_sq = ry * ry, rx * rx
    res = []

    p1=ry_sq-rx_sq*ry+rx_sq/4
    while ry_sq*x<rx_sq*y:
        if f:
            xk, yk = y, x
        else:
            xk, yk = x, y
        res.extend([(xk + xc, yk + yc), (xc - xk, yk + yc),
                    (xk + xc, yc - yk), (xc - xk, yc - yk)])

        if p1 < 0:
            p1 = p1 + 2 * ry_sq * x + 3 * ry_sq
            x = x + 1
        else:
            p1 = p1 + 2 * ry_sq * x - 2 * rx_sq * y + 2 * rx_sq + 3 * ry_sq
            x = x + 1
            y = y - 1

    p2 = ry_sq * (x + 0.5) * (x + 0.5) + rx_sq * \
        (y - 1) * (y - 1) - rx_sq * ry_sq
    while y > 0:
        if f:
            xk, yk = y, x
        else:
            xk, yk = x, y
        res.extend([(xk + xc, yk + yc), (xc - xk, yk + yc),
                    (xk + xc, yc - yk), (xc - xk, yc - yk)])
        if p2 >= 0:
            p2 = p2 - 2 * rx_sq * y + 3 * rx_sq
            y = y - 1
        else:
            p2 = p2 + 2 * ry_sq * x - 2 * rx_sq * y + 2 * ry_sq + 3 * rx_sq
            x = x + 1
            y = y - 1

    if f:
        res.extend([(xc, rx + yc), (xc, yc - rx)])
    else:
        res.extend([(rx + xc, yc), (xc - rx, yc)])
    return res

def Bezier_Point(t, p_list):
    """
    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...])
    :param t: (float)
    :return: (list of int:[x,y])
    """
    new_list = []
    while (len(p_list) > 1):
        for i in range(0, len(p_list) - 1):
            Qx = (1 - t) * p_list[i][0] + t * p_list[i + 1][0]
            Qy = (1 - t) * p_list[i][1] + t * p_list[i + 1][1]
            new_list.append([Qx, Qy])
        p_list = new_list
        new_list = []
    x = int(p_list[0][0])
    y = int(p_list[0][1])
    return x, y


def Basefunction(i, k, u):
    """
    :param i:(int)index of base function
    :param k:(int) degree+1
    :param u:parameter
    :return:the value of base function
    """
    Nik_u = 0.0
    if k == 1:
        if u < i + 1 and u >= i:
            Nik_u = 1.0
        else:
            Nik_u = 0.0
    else:
        Nik_u = ((u - i) / (k - 1)) * Basefunction(i, k - 1, u) + \
            ((i + k - u) / (k - 1)) * Basefunction(i + 1, k - 1, u)

    return Nik_u

def draw_curve(p_list, algorithm):
    """
    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...])
    :param algorithm: (string) 'Bezier' and 'B-spline'
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...])
    """
    res = []
    if algorithm == "Bezier":
        t_step = 0.0005
        t = 0
        while t <= 1:
            res.append(Bezier_Point(t, p_list))
            t = t + t_step
    elif algorithm == "B-spline":
        k = 3
        n = len(p_list) - 1
        u = k
        step = 0.001
        while u <= n + 1:
            p_x = 0.0
            p_y = 0.0
            for i in range(n + 1):
                Nik = Basefunction(i, k + 1, u)
                p_x = p_x + p_list[i][0] * Nik
                p_y = p_y + p_list[i][1] * Nik
            u = u + step
            res.append([int(p_x), int(p_y)])
    return res

def translate(p_list,dx,dy):
    '''
    :param p_list: (list of list of int)
    :param dx:(int)
    :param dy:(int)
    :return:(list of list of int)
    '''
    result=[]
    for i in range(len(p_list)):
        x,y=p_list[i]
        x,y=x+dx,y+dy
        result.append(x,y)
    return result

def rotate(p_list, x, y, r):
    """
    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...])
    :param x: (int)
    :param y: (int)
    :param r: (int)
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...])
    """
    res = []
    r = math.radians(360 + r)
    for x0, y0 in p_list:
        new_x = x + (x0 - x) * math.cos(r) - (y0 - y) * math.sin(r)
        new_y = y + (x0 - x) * math.sin(r) + (y0 - y) * math.cos(r)
        res.append([round(new_x), round(new_y)])
    return res


def scale(p_list, x, y, s):
    """
    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...])
    :param x: (int)
    :param y: (int)
    :param s: (float)
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...])
    """
    res = []
    for x0, y0 in p_list:
        new_x = x + (x0 - x) * s
        new_y = y + (y0 - y) * s
        res.append([round(new_x), round(new_y)])
    return res


def encode(x_min, y_min, x_max, y_max, x, y):
    LEFT   = 0b0001
    RIGHT  = 0b0010
    BOTTOM = 0b0100
    TOP    = 0b1000
    code   = 0
    if x < x_min:
        code += LEFT
    elif x > x_max:
        code += RIGHT
    if y > y_max:
        code += BOTTOM
    elif y < y_min:
        code += TOP
    return code

def clip(p_list, x_min, y_min, x_max, y_max, algorithm):
    """
    :param p_list: (list of list of int: [[x0, y0], [x1, y1]])
    :param x_min:
    :param y_min:
    :param x_max:
    :param y_max:
    :param algorithm: (string) 'Cohen-Sutherland' and 'Liang-Barsky'
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1]])
    """
    if len(p_list) != 2:
        return []
    x0, y0 = p_list[0]
    x1, y1 = p_list[1]
    res = []
    if algorithm == 'Cohen-Sutherland':
        LEFT   = 0b0001
        RIGHT  = 0b0010
        BOTTOM = 0b0100
        TOP    = 0b1000
        code   = [0, 0]
        code[0] = encode(x_min, y_min, x_max, y_max, x0, y0)
        code[1] = encode(x_min, y_min, x_max, y_max, x1, y1)
        res = p_list
        while code[0] | code[1] != 0:
            if code[0] & code[1] != 0:
                return []  # empty
            for i in range(2):
                x, y = res[i]
                if code[i] == 0:
                    continue
                else:
                    if LEFT & code[i] != 0:
                        x = x_min
                        if x0 == x1:
                            y = y0
                        else:
                            y = y0 + (y1 - y0) * (x_min - x0) / (x1 - x0)
                    elif RIGHT & code[i] != 0:
                        x = x_max
                        if x0 == x1:
                            y = y0
                        else:
                            y = y0 + (y1 - y0) * (x_max - x0) / (x1 - x0)
                    elif BOTTOM & code[i] != 0:
                        y = y_max
                        x = x0 + (x1 - x0) * (y_max - y0) / (y1 - y0)
                    elif TOP & code[i] != 0:
                        y = y_min
                        x = x0 + (x1 - x0) * (y_min - y0) / (y1 - y0)
                    code[i] = encode(x_min, y_min, x_max, y_max, x, y)
                    res[i] = [int(x), int(y)]

    elif algorithm == 'Liang-Barsky':
        dx = x1 - x0
        dy = y1 - y0
        p = [-dx, dx, -dy, dy]
        q = [x0 - x_min, x_max - x0, y0 - y_min, y_max - y0]
        u0, u1 = 0, 1
        for k in range(4):
            if p[k] == 0:
                if q[k] < 0:
                    return []
            else:
                u = q[k] / p[k]
                if p[k] < 0:
                    u0 = max(u0, u)
                else:
                    u1 = min(u1, u)
        if u0 > u1:
            return []
        x_0 = round(x0 + u0 * dx)
        y_0 = round(y0 + u0 * dy)
        x_1 = round(x0 + u1 * dx)
        y_1 = round(y0 + u1 * dy)
        res = [[x_0, y_0], [x_1, y_1]]

    return res
