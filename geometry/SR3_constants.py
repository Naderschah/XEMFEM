import math

liquid_level = 0.004
### ----------------------  Helper

def transform_pts(pts, dx=0.0, dy=0.0, mirror_y=False):
    """
    pts items:
      ["line", x, y]
      ["arc",  x, y, cx, cy]
    Applies:
      x := x + dx
      y := (+/-) y + dy   (mirror_y controls sign)
    Also transforms arc centers consistently.
    """
    out = []
    for item in pts:
        tag = item[0]

        if tag == "line":
            x, y = item[1], item[2]
            y = (-y if mirror_y else y)
            out.append(["line", x + dx, y + dy])

        elif tag == "arc":
            if (len(item) != 5) and (len(item) != 6):
                raise ValueError(f"Bad arc item length: {item}")
            if len(item) == 5: direction = True
            if len(item) == 6: direction = item[5]
            x, y, cx, cy = item[1], item[2], item[3], item[4]
            y  = (-y  if mirror_y else y)
            cy = (-cy if mirror_y else cy)
            out.append(["arc", x + dx, y + dy, cx + dx, cy + dy, direction])

        else:
            raise ValueError(f"Unknown pts tag: {tag}")

    return out


### ----------------------  PTFE Elements
def PTFEWall(in_dict):
    H = 1.5008
    part_dict = {
        'Width': 0.003,
        'Height': H,
        'VerticalPosition': -0.0008 - H,
        'RadialPosition': 0.664
    }
    in_dict['PTFEWall'] = part_dict
    return in_dict

def CopperRingInsulatingFrame(in_dict):
    H = 0.005
    part_dict = {
        'Width': 0.028,
        'Height': H,
        'VerticalPosition': -0.02-H,
        'RadialPosition': 0.679
    }
    in_dict['CopperRingInsulatingFrame'] = part_dict
    return in_dict

def BottomScreenInsulatingFrame(in_dict):
    part_dict = {
        'Width': 0.002,
        'Height': 0.005,
        'FilletRadius': 0.001,
        'VerticalPosition': 0.0407,
        'RadialPosition': -1.558
    }
    in_dict['BottomScreenInsulatingFrame'] = part_dict
    return in_dict

def TopScreenInsulatingFrame(in_dict):
    # ---- Inherent Geometry
    A = 0.0165
    B = 0.0005
    C = 0.0035
    D = 0.0335
    E = 0.0149
    F = 0.0027
    G = 0.0332
    H = 0.0093
    I = 0.0222
    J = 0.002

    # ---- Global Alignment
    VerticalPosition = 0.0407
    RadialPosition = 0.6643

    # ---- X stations
    x0 = RadialPosition
    x1 = x0 + F
    x2 = x1 + D
    x3 = x2 - J
    x4 = x3 - I

    # ---- Y stations
    y0 = VerticalPosition
    y1 = y0 + G
    y2 = y1 - E
    y3 = y2 - C
    y4 = y3 - B
    y5 = y4 - A

    part_dict = {
        'A': A,
        'B': B,
        'C': C,
        'D': D,
        'E': E,
        'F': F,
        'G': G,
        'H': H,
        'I': I,
        'J': J,
        'VerticalPosition': VerticalPosition,
        'RadialPosition': RadialPosition,
        'pts': [
            [x0, y1],
            [x1, y1],
            [x1, y2],
            [x2, y2],
            [x2, y5],
            [x3, y5],
            [x3, y4],
            [x4, y4],
            [x4, y3],
            [x1, y3],
            [x0, y3],
        ]
    }

    in_dict['TopScreenInsulatingFrame'] = part_dict
    return in_dict

def CathodeInsulatingFrame(in_dict):
    # ---- Inherent Geometry
    A = 0.0265
    B = 0.019
    C = 0.013
    D = 0.0055
    E = 0.0045

    F = 0.00114
    G = 0.003
    H = 0.0275
    I = 0.002
    J = 0.002
    K = 0.018

    VerticalOffset = 0.00036

    # ---- Global Alignment
    VerticalPosition = -1.5028
    RadialPosition = 0.668 # FIXME Is it 0.668 - 0.66914 most left corner in old script

    # ---- X stations
    x0 = RadialPosition + F
    x1 = x0 + G
    x2 = x1 + H
    x3 = x2 - I
    x4 = x3 - J
    x5 = x4 - K

    # ---- Y stations (pre-shrink)
    y0 = VerticalPosition + A
    y1 = y0 - B
    y2 = y1 - C
    y3 = y2 + D
    y4 = y3 + E

    # ---- Final offset station
    ys4o = y4 + VerticalOffset

    part_dict = {
        'pts': [
            [x0, y0],
            [x1, y0],
            [x1, y1],
            [x2, y1],
            [x2, y2],
            [x3, y2],
            [x3, y3],
            [x4, y3],
            [x4, ys4o],
            [x5, ys4o],
            [x5, y3],
            [x0, y3],
        ]
    }

    in_dict['CathodeInsulatingFrame'] = part_dict
    return in_dict

def BottomStackInsulatingFrame(in_dict):
    """
    This one is not fun, we have a rectangle with holes cut in it
    """
    # --- The thing in front of the bottom electrodes
    VerticalPosition = -1.5035
    GrooveVerticalPos = 0.00178
    W = 0.0042
    H = 0.05378
    x0 = 0.664
    y0 = VerticalPosition - 0.05378
    y_top = y0 + H
    r = 0.001 / 2
    pitch = 0.002
    n_rows = int(H / pitch)

    # groove center y positions (top row first)
    yc = [VerticalPosition - GrooveVerticalPos - i * pitch for i in range(n_rows)]

    x_left = x0
    x_right = x0 + W
    yc_first = VerticalPosition - GrooveVerticalPos
    yc_final = VerticalPosition - GrooveVerticalPos - (n_rows - 1) * pitch

    # ---- pts (verbatim logic)
    pts = []
    for i in range(n_rows):
        yc_i = yc_first - i * pitch
        # arc start vertex (connects by arc to the next vertex)
        pts.append(["arc", x_left, yc_i + r, x_left, yc_i])
        # arc end vertex (connects by line to next arc start, or to bottom-left corner if last)
        pts.append(["line", x_left, yc_i - r])

    # -------------------------
    # BOTTOM (left -> right)
    # -------------------------
    # After the last left groove end: line to bottom-left corner, then to bottom-right corner,
    # then to the start of the bottom-right groove arc (arcs_right[0].start).
    y_bot = (-1.5035 - 0.05378)
    pts.append(["line", x_left,  y_bot])           # bottom-left corner
    pts.append(["line", x_right, y_bot])           # bottom-right corner

    # -------------------------
    # RIGHT SIDE (bottom -> top)
    # -------------------------
    # bottommost right arc start vertex (connects by arc to its end)
    pts.append(["arc", x_right, yc_final - r, x_right, yc_final])
    for k in range(n_rows - 1, 0, -1):
        # current arc center
        yc_k = yc_first - k * pitch          # this equals yc[k]
        yc_km1 = yc_first - (k - 1) * pitch  # next one up
        # end of current arc (connects by line to next arc start)
        pts.append(["line", x_right, yc_k + r])
        # start of next arc up (connects by arc to its end)
        pts.append(["arc", x_right, yc_km1 - r, x_right, yc_km1])

    pts.append(["line", x_right, yc_first + r])
    pts.append(["line", x_right, y_top])          # top-right corner
    pts.append(["line", x_left,  y_top])

    part_dict = {
        'pts': pts,
    }

    in_dict['BottomStackInsulatingFrame'] = part_dict
    return in_dict

def PMTReflectorShared():
    BottomVerticalPosition = -1.56
    TopVerticalPosition = 0.0739
    x0 = 0.648

    P = 0.081
    a = 0.5 * P
    P_fix = 0.081

    xh = 0.0785 / 2.0
    xp = 0.069 / 2.0
    xc = xp - 0.002

    yT = 0.0
    yC = -0.002
    yD = -0.0029
    yB = -0.008

    offset_x = 0.0785

    pts = [
        [a - xh, yB],
        [a - xh, yD],
        [a - xc, yD],
        [a - xc, yC],
        [a - xp, yT],
        [-a + xp, yT],
        [-a + xc, yC],
        [-a + xc, yD],
        [-a + xh, yD],
        [-a + xh, yB],
    ]

    return {
        'BottomVerticalPosition': BottomVerticalPosition,
        'TopVerticalPosition': TopVerticalPosition,
        'x0': x0,
        'P': P,
        'a': a,
        'xh': xh,
        'xp': xp,
        'xc': xc,
        'yT': yT,
        'yC': yC,
        'yD': yD,
        'yB': yB,
        'pts': pts,
        'offset_x': offset_x,
        "P_fix": P_fix
    }

def BottomPMTReflectorRightBit(in_dict):
    s = PMTReflectorShared()

    BottomVerticalPosition = s['BottomVerticalPosition']
    x0 = s['x0']
    xp = s['xp']
    xh = s['xh']

    yT = BottomVerticalPosition
    xc = xp - 0.002


    yC = yT - 0.002
    yD = yT - 0.0029
    yB = yT - 0.008

    xR  = 1.395 / 2.0
    xN  = 0.69025
    xNW = xN + 0.002
    yN  = yT - 0.003

    pts = [
        [x0 + xp, yT],
        [x0 + xc, yC],
        [x0 + xc, yD],
        [x0 + xh, yD],
        [x0 + xh, yB],
        [xR,      yB],
        [xR,      yT],
        [xNW,     yT],
        [xNW,     yN],
        [xN,      yN],
        [xN,      yT],
    ]

    part_dict = {
        'pts': pts,
    }

    in_dict['BottomPMTReflectorRightBit'] = part_dict
    return in_dict

def BottomPMTReflectors(in_dict):
    s = PMTReflectorShared()

    BottomVerticalPosition = s['BottomVerticalPosition']
    P = s['P']
    P_fix = s['P_fix']
    base_pts = s['pts']

    yT = BottomVerticalPosition
    # Only this makes the midpoint intercept at x = 0
    x0 = s['x0'] - P/2

    pts = [[i[0] + x0, i[1] + yT] for i in base_pts]

    part_dict = {
        'pts': pts,
        'HorizontalPitch': -P,
        'Number': 8,
    }

    in_dict['BottomPMTReflectors'] = part_dict
    return in_dict

def TopPMTReflectors(in_dict):
    s = PMTReflectorShared()

    P = s['P']
    P_fix = s['P_fix']
    base_pts = s['pts']

    TopVerticalPosition = s['TopVerticalPosition']
    yT = TopVerticalPosition
    x0 = s['x0'] - P/2

    # NOTE: this follows your class verbatim: i[0] + x - P (uses x, not x0)
    pts = [[i[0] + x0, -i[1] + yT] for i in base_pts]

    part_dict = {
        'pts': pts,
        'HorizontalPitch': -P,
        'Number': 8,
    }

    in_dict['TopPMTReflectors'] = part_dict
    return in_dict

def TopPMTReflectorRightBit(in_dict):
    s = PMTReflectorShared()

    # shared
    TopVerticalPosition = s["TopVerticalPosition"]
    x0 = s["x0"]
    xp = s["xp"]
    xh = s["xh"]
    yC_shared = s["yC"]
    yD_shared = s["yD"]
    yB_shared = s["yB"]

    # class verbatim
    P = 0.081

    # half-widths and chamfer inset
    xc = xp + yC_shared  # Note yC gets redefined below

    # vertical references (top row uses + height/depth above VerticalPosition)
    yT = TopVerticalPosition
    yC = yT - yC_shared
    yD = yT - yD_shared
    yB = yT - yB_shared

    # x references
    xR = 1.412 / 2.0

    pts = [
        [x0 + xp, yT],  # hole opening at top surface
        [x0 + xc, yC],  # chamfer up/out (top reflector geometry is "above" y0)
        [x0 + xc, yD],  # vertical to hole depth
        [x0 + xh, yD],  # step out to housing at depth
        [x0 + xh, yB],  # outer housing up to reflector height
        [xR,      yB],  # outer reflector boundary at height
        [xR,      yT],  # outer reflector boundary back down to top surface
    ]
    part_dict = {
        "pts": pts,
    }
    in_dict["TopPMTReflectorRightBit"] = part_dict
    return in_dict


def PMTReflectors(in_dict):
    sub_dict = {}

    sub_dict = BottomPMTReflectors(sub_dict)
    sub_dict = BottomPMTReflectorRightBit(sub_dict)
    sub_dict = TopPMTReflectors(sub_dict)
    sub_dict = TopPMTReflectorRightBit(sub_dict)

    part_dict = {
        'hull': True,
        'sub_sketches': sub_dict,
    }

    in_dict['PMTReflectors'] = part_dict
    return in_dict

def UnderCathodePTFEPin(in_dict):
    yplus = 0.005
    y0 = -1.558 - yplus
    W = 0.002
    x0 = 0.69025

    pts = [
        [x0, y0],
        [x0 + W, y0],
        [x0 + W, y0 + yplus],
        [x0, y0 + yplus],
    ]

    part_dict = {
        'pts': pts,
    }

    in_dict['UnderCathodePTFEPin'] = part_dict
    return in_dict

def AnodeInsulatingFrame(in_dict):
    # ---- Inherent Geometry
    A = 0.0308
    B = 0.0005
    C = 0.012
    D = 0.0202
    E = 0.004
    F = 0.002
    G = 0.0315
    H = 0.002
    I = 0.022
    J = 0.0005
    K = 0.0027
    N = 0.0148
    L = 0.0255

    # ---- Global Alignment
    VerticalPosition = 0.0087
    RadialPosition = 0.6643

    # ---- X stations
    x0 = RadialPosition
    x1 = x0 + K
    x2 = x0 + C
    x3 = x2 + D
    x4 = x3 + E
    x5 = x4 - H
    x6 = x5 - I

    # ---- Y stations
    y0 = VerticalPosition
    y1 = y0 + A
    y2 = y1 + N
    y3 = y1 + B
    y4 = y3 - F
    y5 = y4 - G
    y6 = y5 + L
    y7 = y6 + J

    pts = [
        [x0, y0],
        [x0, y1],
        [x0, y2],
        [x1, y2],
        [x1, y1],
        [x2, y1],
        [x2, y3],
        [x3, y3],
        [x3, y4],
        [x4, y4],
        [x4, y5],
        [x5, y5],
        [x5, y6],
        [x6, y6],
        [x6, y7],
        [x1, y7],
        [x1, y0],
    ]

    part_dict = {
        'pts': pts,
    }

    in_dict['AnodeInsulatingFrame'] = part_dict
    return in_dict

def GateInsulatingFrame(in_dict):
    # ---- Inherent Geometry
    A = 0.026
    B = 0.002
    C = 0.004
    D = 0.0202
    E = 0.012
    F = 0.0005
    G = 0.007
    H = 0.0005
    I = 0.0222
    J = 0.02
    K = 0.002

    # ---- Global Alignment
    VerticalPosition = 0.0005
    RadialPosition = 0.6643

    # ---- X stations
    x0 = RadialPosition
    x1 = x0 + E
    x2 = x1 + D
    x3 = x2 + C
    x4 = x3 - K

    # ---- Y stations
    y0 = VerticalPosition
    y1 = y0 + G
    y2 = y1 + F
    y3 = y2 - B
    y4 = y3 - A
    y5 = y4 + J
    y6 = y0 - H

    pts = [
        [x0, y0],
        [x0, y1],
        [x1, y1],
        [x1, y2],
        [x2, y2],
        [x2, y3],
        [x3, y3],
        [x3, y4],
        [x4, y4],
        [x4, y5],
        [x1, y6],
        [x1, y0],
    ]

    part_dict = {
        'pts': pts,
    }

    in_dict['GateInsulatingFrame'] = part_dict
    return in_dict


### -----------------------  PMTs ----------------------------------------------


def PMTsShared():
    # Defines the PMT as symmetric about the x axis with the window touching the
    # ---------- Variables used as helpers
    QD = 0.076      # Quartz Diameter
    QH = 0.0135     # Quartz constants
    BD = 0.0533     # Body Diameter
    BH = 0.1005     # Body Height
    BaD = 0.0318    # PMT base diameter
    BaH = 0.00955   # PMT base height

    # ---------- Variables for Geometry creation
    # Quartz parts
    v0x = -QD / 2
    v0y = 0
    v1x = +QD / 2
    v2y = -QH

    # Bell Body Intersection
    v3x = BD / 2
    v3y = -QH - math.sqrt((QD / 2) ** 2 - (BD / 2) ** 2)
    v4x = +BD / 2
    v4y = -QH - BH

    # Base
    v5x = +BaD / 2
    v6x = +BaD / 2
    v6y = -(QH + BH + BaH)
    v7x = -BaD / 2
    v8y = -QH - BH

    # Body Left
    v9x = -BD / 2

    # Bell Body intersection
    v10x = -BD / 2

    # Quartz bottom left
    v11x = -QD / 2

    pts = [
        ["line", v0x,  v0y],                 # 0 -> 1  quartz face (top)
        ["line", v1x,  v0y],                 # 1 -> 2  quartz right wall
        ["arc",  v1x,  v2y, 0, -QH, True],         # 2 -> 3  right bell arc
        ["line", v3x,  v3y],                 # 3 -> 4  body right wall
        ["line", v4x,  v4y],                 # 4 -> 5  body to base shoulder
        ["line", v5x,  v4y],                 # 5 -> 6  base right wall
        ["line", v6x,  v6y],                 # 6 -> 7  base bottom
        ["line", v7x,  v6y],                 # 7 -> 8  base left wall
        ["line", v7x,  v8y],                 # 8 -> 9  base to body shoulder
        ["line", v9x,  v8y],                 # 9 -> 10 body left wall
        ["arc",  v10x, v3y, 0, -QH, True],         # 10 -> 11 left bell arc
        ["line", v11x, v2y],                 # 11 -> 0  quartz left wall (closes loop)
    ]

    return {
        'QD': QD, 'QH': QH, 'BD': BD, 'BH': BH, 'BaD': BaD, 'BaH': BaH,
        'pts': pts, 
        'pitch': 0.648/8
    }

def InnerPMTBottom(in_dict):
    s = PMTsShared()
    pts = s['pts']

    VerticalPosition = -1.56 - 0.0029
    RadialPosition = 0

    # Shift to location we care for bottom x = 0 PMT
    pts = transform_pts(pts, dx=0.0, dy=VerticalPosition, mirror_y=False)

    # pts[0] is the negative x extent of the outer portion of the quartz
    pts[0][1] = 0.0
    # pts[7] is the most negative x extent of the base portion of the PMT
    pts[7][1] = 0.0

    # Everything else can be trashed
    pts = pts[0:8]

    part_dict = {
        'pts': pts,
        'VerticalPosition': VerticalPosition,
        'RadialPosition': RadialPosition,
    }

    in_dict['InnerPMTBottom'] = part_dict
    return in_dict

def BottomPMTs(in_dict):
    s = PMTsShared()
    pts = s['pts']
    P = s['pitch']

    VerticalPosition = -1.56 - 0.0029
    RadialPosition = 0.648

    # Shift to location we care for bottom PMTs
    pts = transform_pts(pts, dx=RadialPosition, dy=VerticalPosition, mirror_y=False)

    sub_dict = {}
    sub_dict = InnerPMTBottom(sub_dict)

    part_dict = {
        'pts': pts,
        'HorizontalPitch': -P,
        'Number': 8,
        'sub_sketches': sub_dict,
        "RadialPosition":RadialPosition,
        "VerticalPosition":VerticalPosition
    }

    in_dict['BottomPMTs'] = part_dict
    return in_dict

def InnerPMTTop(in_dict):
    s = PMTsShared()
    pts = s['pts']

    VerticalPosition = 0.0739 + 0.0029
    RadialPosition = 0

    # Shift to location we care for top x = 0 PMT (mirrored in y)
    pts = transform_pts(pts, dx=0.0, dy=VerticalPosition, mirror_y=True)

    # pts[0] is the negative x extent of the outer portion of the quartz
    pts[0][1] = 0.0
    # pts[7] is the most negative x extent of the base portion of the PMT
    pts[7][1] = 0.0
    # Need to change direction of arc for top PMTs
    pts[2][-1] = False

    # Everything else can be trashed
    pts = pts[0:8]

    part_dict = {
        'pts': pts,
        'VerticalPosition': VerticalPosition,
        'RadialPosition': RadialPosition,
    }

    in_dict['InnerPMTTop'] = part_dict
    return in_dict

def TopPMTs(in_dict):
    s = PMTsShared()
    pts = s['pts']
    P = s['pitch']

    VerticalPosition = 0.0739 + 0.0029
    RadialPosition = 0.648

    # Shift to location we care for top PMTs (mirrored in y)
    pts = transform_pts(pts, dx=RadialPosition, dy=VerticalPosition, mirror_y=True)

    # Need to change direction of arc for top PMTs
    pts[2][-1] = False
    pts[-2][-1] = False


    sub_dict = {}
    sub_dict = InnerPMTTop(sub_dict)

    part_dict = {
        'pts': pts,
        'HorizontalPitch': -P,
        'Number': 8,
        'sub_sketches': sub_dict,
    }

    in_dict['TopPMTs'] = part_dict
    return in_dict

def AllPMTs(in_dict):
    sub_dict = {}

    sub_dict = TopPMTs(sub_dict)
    sub_dict = BottomPMTs(sub_dict)

    part_dict = {
        'hull': True,
        'sub_sketches': sub_dict,
    }

    in_dict['AllPMTs'] = part_dict
    return in_dict

# -------------------------  Field Shaping ------------------------------------------
def FieldShapingRings_5(in_dict):
    Radius = 0.001
    VerticalPosition = -0.023
    RadialPosition = 0.668
    VerticalPitch = 0.022 / 2
    Number = 5

    part_dict = {
        'Radius': Radius,
        'VerticalPosition': VerticalPosition,
        'RadialPosition': RadialPosition,
        'VerticalPitch': -VerticalPitch,
        'Number': Number,
    }

    in_dict['FieldShapingRings_5'] = part_dict
    return in_dict


def FieldShapingRings_65(in_dict, prev_vertical_position):
    Radius = 0.001
    RadialPosition = 0.668
    VerticalPitch = 0.022

    VerticalPosition = prev_vertical_position - 4 * (VerticalPitch / 2)
    Number = 65

    part_dict = {
        'Radius': Radius,
        'VerticalPosition': VerticalPosition,
        'RadialPosition': RadialPosition,
        'VerticalPitch': -VerticalPitch,
        'Number': Number,
    }

    in_dict['FieldShapingRings_65'] = part_dict
    return in_dict


def FieldShapingRings_2(in_dict, prev_vertical_position):
    Radius = 0.001
    RadialPosition = 0.668
    VerticalPitch = 0.022 / 2

    VerticalPosition = prev_vertical_position - 64.5 * (2 * VerticalPitch)
    Number = 2

    part_dict = {
        'Radius': Radius,
        'VerticalPosition': VerticalPosition,
        'RadialPosition': RadialPosition,
        'VerticalPitch': -VerticalPitch,
        'Number': Number,
    }

    in_dict['FieldShapingRings_2'] = part_dict
    return in_dict


def FieldShapingRings(in_dict):
    sub_dict = {}

    sub_dict = FieldShapingRings_5(sub_dict)
    vp_5 = sub_dict['FieldShapingRings_5']['VerticalPosition']

    sub_dict = FieldShapingRings_65(sub_dict, vp_5)
    vp_65 = sub_dict['FieldShapingRings_65']['VerticalPosition']

    sub_dict = FieldShapingRings_2(sub_dict, vp_65)

    part_dict = {
        'hull': True,
        'sub_sketches': sub_dict,
    }

    in_dict['FieldShapingRings'] = part_dict
    return in_dict


def FieldShapingGuard(in_dict):
    H = 0.015
    W = 0.005
    y0 = -0.078125 * in_dict['shrinkage_factor'] - H /2
    x0 = 0.677687

    r = W / 2.00

    pts = [
        # 0) bottom semicircle: left -> right
        ["arc",  x0,      y0 + r,      x0 + W/2.0,  y0 + r,      False],
        # 1) right vertical: bottom -> top
        ["line", x0 + W,  y0 + r],
        # 2) top semicircle: right -> left
        ["arc",  x0 + W,  y0 + H - r,  x0 + W/2.0,  y0 + H - r,  False],
        # 3) left vertical: top -> bottom (closes)
        ["line", x0,      y0 + H - r],
    ]

    part_dict = {
        #"Width": W,
        #"Height": H,
        #"VerticalPosition": y,
        #"RadialPosition": x,
        'pts':pts,
        #"FilletRadius": Rf,
        #"FilletIndeces": [0,1,2,3],
        # omit FilletIndeces to mean "all edges" (per your template comment)
        "VerticalPitch": -0.022,
        "Number": 64,
    }

    in_dict["FieldShapingGuard"] = part_dict
    return in_dict

# -----------------------------  Electrodes ----------------------------------------------
def GateWires(in_dict):
    VerticalPosition = 0
    Radius = 0.000216
    HorizontalPitch = 0.005
    RadialPosition = HorizontalPitch / 4
    Number = 133

    part_dict = {
        'VerticalPosition': VerticalPosition,
        'Radius': Radius,
        'HorizontalPitch': HorizontalPitch,
        'RadialPosition': RadialPosition,
        'Number': Number,
    }

    in_dict['GateWires'] = part_dict
    return in_dict
def GateElectrode(in_dict):
    Height = 0.02
    Width = 0.031
    CutOutHeight = 0.011
    CutOutWidth = 0.01
    FilletRadius = 0.0016

    VerticalPosition = 0
    RadialPosition = 0.667

    sub_dict = {}
    sub_dict = GateWires(sub_dict)

    x0 = RadialPosition
    y0 = VerticalPosition
    H  = Height
    W  = Width
    ch = CutOutHeight
    cw = CutOutWidth

    pts = [
        ["line", x0,      y0 - H + ch],
        ["line", x0,      y0],
        ["line", x0 + W,  y0],
        ["line", x0 + W,  y0 - H],
        ["line", x0 + cw, y0 - H],
        ["line", x0 + cw, y0 - H + ch],
    ]
    part_dict = {
        "pts": pts,
        "FilletRadius": FilletRadius,   
        "FilletIndeces": [0, 1, 2, 3, 4], 
        "sub_sketches": sub_dict,
    }

    in_dict["GateElectrode"] = part_dict
    return in_dict

def CathodeWires(in_dict):
    VerticalPosition = -1.5028

    Radius = 0.0003
    HorizontalPitch = -0.0075
    RadialPosition = 0.6675 - 3/4 * HorizontalPitch
    Number = 89

    part_dict = {
        'VerticalPosition': VerticalPosition,
        'Radius': Radius,
        'HorizontalPitch': HorizontalPitch,
        'RadialPosition': RadialPosition + HorizontalPitch,
        'Number': Number,
    }

    in_dict['CathodeWires'] = part_dict
    return in_dict
def CathodeElectrode(in_dict):
    x0 = 0.6735
    y0 = (-1.5028 - 0.02)
    W = 0.024
    H = 0.02

    rR = 0.01    # for A, J
    rD = 0.0005    # for B, I
    rF = 0.001                       # for F
    rG = 0.002

    pts =  [
        # --- A fillet (bottom-left) is between bottom edge and left edge ---
        ["line", x0,          y0 + rR],                 # A_out: up left wall
        #["line", x0,          y0 + H - rD],             # B_in (before top-left fillet)

        # --- B fillet (top-left) between left wall and top edge ---
        ["arc",  x0,          y0 + H - rD,  x0 + rD,    y0 + H - rD],   # B_in -> B_out
        ["line", x0 + rD,     y0 + H],                  # B_out: along top edge to C

        # --- C, D, E are sharp corners (no fillets) ---
        ["line", x0 + 0.0045,                 y0 + H],                  # C
        ["line", x0 + 0.0045,                 y0 + H + 0.003],          # D
        ["line", x0 + 0.0045 + 0.002,         y0 + H + 0.003],          # E

        # --- F fillet (radius rF) between vertical up (E->F) and horizontal right (F->G) ---
        #["line", x0 + 0.0065,                 y0 + H + 0.0048 - rF],    # F_in
        ["arc",  x0 + 0.0065,                 y0 + H + 0.0048 - rF,
                x0 + 0.0065 + rF,            y0 + H + 0.0048 - rF],    # F_in -> F_out
        ["line", x0 + 0.0065 + rF,            y0 + H + 0.0048],         # F_out

        # --- G fillet (radius rG) between horizontal right (F->G) and vertical down (G->H) ---
        #["line", x0 + 0.0195 - rG,            y0 + H + 0.0048],         # G_in
        ["arc",  x0 + 0.0195 - rG,            y0 + H + 0.0048,
                x0 + 0.0195 - rG,            y0 + H + 0.0048 - rG],    # G_in -> G_out
        ["line", x0 + 0.0195,                 y0 + H + 0.0048 - rG],    # G_out: down toward H

        # --- H and I corner region (H is sharp; I is filleted) ---
        ["line", x0 + 0.0195,                 y0 + H],                  # H
        #["line", x0 + W - rD,                 y0 + H],                  # I_in

        # --- I fillet (top-right) between top edge and right wall ---
        ["arc",  x0 + W - rD,                 y0 + H,  x0 + W - rD,      y0 + H - rD],   # I_in -> I_out
        ["line", x0 + W,                      y0 + H - rD],             # I_out: down right wall

        # --- J fillet (bottom-right) between right wall and bottom edge ---
        #["line", x0 + W,                      y0 + rR],                 # J_in
        ["arc",  x0 + W,                      y0 + rR,  x0 + W - rR,     y0 + rR],       # J_in -> J_out
        ["line", x0 + W - rR,                 y0],                      # J_out: along bottom edge

        # --- close back to start (A_out) ---
        #["line", x0 + rR,                     y0],                      # A_in
        ["arc",  x0 + rR,                     y0,      x0 + rR,          y0 + rR],       # A_in -> A_out (closes)
    ]


    sub_dict = {}
    sub_dict = CathodeWires(sub_dict)

    part_dict = {
        'pts': pts,
        'sub_sketches': sub_dict,
    }

    in_dict['CathodeElectrode'] = part_dict
    return in_dict

def AnodeWires(in_dict):
    Radius = 0.000216 / 2
    VerticalPosition = 0.008
    HorizontalPitch = 0.005
    RadialPosition = HorizontalPitch / 4
    Number = 133

    part_dict = {
        'Radius': Radius,
        'VerticalPosition': VerticalPosition,
        'HorizontalPitch': HorizontalPitch,
        'RadialPosition': RadialPosition,
        'Number': Number,
    }

    in_dict['AnodeWires'] = part_dict
    return in_dict
def AnodeElectrode(in_dict):
    Height = 0.024
    Width = 0.031
    FilletRadius = 0.0016

    VerticalPosition = 0.008
    RadialPosition = 0.667

    sub_dict = {}
    sub_dict = AnodeWires(sub_dict)

    part_dict = {
        'Height': Height,
        'Width': Width,
        'FilletRadius': FilletRadius,
        'VerticalPosition': VerticalPosition,
        'RadialPosition': RadialPosition,
        'sub_sketches': sub_dict,
    }

    in_dict['AnodeElectrode'] = part_dict
    return in_dict

# ------------------- Screening Electrodes 
def TopScreeningElectrode(in_dict):
    Height = 0.015
    Width = 0.031
    FilletRadius = 0.0016

    VerticalPosition = 0.055
    RadialPosition = 0.667

    part_dict = {
        'Height': Height,
        'Width': Width,
        'FilletRadius': FilletRadius,
        'VerticalPosition': VerticalPosition-Height,
        'RadialPosition': RadialPosition,
    }

    in_dict['TopScreeningElectrode'] = part_dict
    return in_dict

def BottomScreeningElectrode(in_dict):
    # ---- Inherent Geometry
    H = 0.015
    W = 0.025
    rD = 0.0005
    rR = 0.0075

    # ---- Global Alignment
    y0 = -1.558
    x0 = 0.6725

    pts = [
        # left wall up to start of top-left fillet
        ["line", x0,           y0 + rD],                      # -> (x0, y0+H-rR)
        # top-left fillet
        ["arc",  x0,           y0 + H - rR,  x0 + rR,  y0 + H - rR, True],  # -> (x0+rR, y0+H)
        # top edge to start of top-right fillet
        ["line", x0 + rR,      y0 + H],                       # -> (x0+W-rR, y0+H)
        # top-right fillet
        ["arc",  x0 + W - rR,  y0 + H,      x0 + W - rR, y0 + H - rR, True], # -> (x0+W, y0+H-rR)
        # right wall down to start of bottom-right fillet
        ["line", x0 + W,       y0 + H - rR],                  # -> (x0+W, y0+rD)
        # bottom-right fillet
        ["arc",  x0 + W,       y0 + rD,     x0 + W - rD, y0 + rD, True],      # -> (x0+W-rD, y0)
        # bottom edge to start of bottom-left fillet
        ["line", x0 + W - rD,  y0],                           # -> (x0+rD, y0)
        # bottom-left fillet (closes to first point)
        ["arc",  x0 + rD,      y0,          x0 + rD,     y0 + rD, True],      # -> (x0, y0+rD)
    ]

    part_dict = {
        'pts': pts,
    }

    in_dict['BottomScreeningElectrode'] = part_dict
    return in_dict

def CopperRing(in_dict):
    Height = 0.01
    Width = 0.021
    FilletRadius = 0.0015

    VerticalPosition = -0.025
    RadialPosition = 0.679

    part_dict = {
        'Height': Height,
        'Width': Width,
        'FilletRadius': FilletRadius,
        'VerticalPosition': VerticalPosition-Height,
        'RadialPosition': RadialPosition,
    }

    in_dict['CopperRing'] = part_dict
    return in_dict


def Bell(in_dict):
    x0 = 0.0
    y_top    = 0.239 
    D_outer  = 0.713
    H_bell   = 0.264
    t_minor  = 0.004
    t_major  = 0.005
    y_thick  = 0.159

    pts = [
        [x0,                y_top],                 # 0 top-left
        [D_outer,           y_top],                 # 1 top-right
        [D_outer,           y_top - H_bell],        # 2 right-bottom thin
        [D_outer - t_minor, y_top - H_bell],        # 3 inner at thin
        [D_outer - t_minor, y_thick],               # 4 inner up to thick transition
        [D_outer - t_major, y_thick],               # 5 inner at thick
        [D_outer - t_major, y_top - t_major],       # 6 inner top of thick section
        [x0,                y_top - t_major],       # 7 left inner
    ]

    in_dict['Bell'] = {
        'pts': pts,
    }
    return in_dict

# -------------------- L & GXe ----------------------

def LXe(in_dict):
    y1 = -1.60824
    y2 = 0.004
    R  = 1.46 / 2.0

    R_upper = 0.22
    R_lower = 1.2

    term_radius_diff = R_lower - R_upper
    dx = R - R_upper
    alpha = math.asin(dx / term_radius_diff)

    dy = math.sqrt((R_lower - R_upper) ** 2 - (R - R_upper) ** 2)
    cy11 = y1 + dy

    # --- upper rounding circle (arc10) ---
    cx10 = R - R_upper
    cy10 = y1

    pts = [
        ["line", 0.0, y1],                                   # axis up
        ["line", 0.0, y2],                                   # top out
        ["line", R,   y2],                                   # right wall down

        # upper rounding arc
        ["arc",  R,   y1,  cx10, cy10],

        # lower rounding arc
        [
            "arc",
            R_lower * math.cos(-math.pi / 2.0 + alpha),
            cy11 + R_lower * math.sin(-math.pi / 2.0 + alpha),
            0.0,
            cy11,
        ],

        # closing line
        [
            "line",
            R_lower * math.cos(-math.pi / 2.0),
            cy11 + R_lower * math.sin(-math.pi / 2.0),
        ],
    ]

    part_dict = {
        'pts': pts,
    }

    in_dict['LXe'] = part_dict
    return in_dict

def GXe(in_dict):
    R = 1.46 / 2.0

    y1 = liquid_level
    y2 = 1.8671 + -1.60824

    R_lower_top = 0.22
    R_upper_top = 1.2

    dx = R - R_lower_top
    deltaR = R_upper_top - R_lower_top

    alpha = math.asin(dx / deltaR)
    dy = math.sqrt(deltaR ** 2 - dx ** 2)

    # c13 (small arc) center
    cx13 = 0.0
    cy13 = y2 - dy
    r13 = R_upper_top

    theta_s_small = math.pi / 2.0 - alpha
    theta_e_small = math.pi / 2.0

    # common point with c12
    x_common = cx13 + r13 * math.cos(theta_s_small)
    y_common = cy13 + r13 * math.sin(theta_s_small)

    # axis-side point
    x_axis = cx13 + r13 * math.cos(theta_e_small)   # = 0
    y_axis = cy13 + r13 * math.sin(theta_e_small)

    # c12 (big arc) center
    cx12 = R - R_lower_top
    cy12 = y2

    pts = [
        ["line", 0.0, y1],   # (0,y1) -> (0,y2) axis up
        ["line", 0.0, y2],   # (0,y2) -> (0,y_axis) axis-cap line start is same point; we keep contour order below

        # Instead of duplicating (0,y2), we go along the liquid-level "bottom flat" first, matching your sketch wires:
        ["line", 0.0, y1],   # (0,y1) -> (R,y1) bottom flat at liquid level
        ["line", R,   y1],   # (R,y1) -> (R,y2) right wall up

        # big arc on c12: (R,y2) -> (x_common,y_common) about (cx12,cy12)
        ["arc",  R,   y2,  cx12, cy12, False],

        # small arc on c13: (x_common,y_common) -> (x_axis,y_axis) about (cx13,cy13)
        ["arc",  x_common, y_common,  cx13, cy13, False],

        # axis-cap line: (x_axis,y_axis) -> (0,y2)
        ["line", x_axis, y_axis],

        # axis down: (0,y2) -> (0,y1) closes
        ["line", 0.0, y2],
    ]

    part_dict = {
        'pts': pts,
    }

    in_dict['GXe'] = part_dict
    return in_dict


# ---------------- Exported Function ----------------
def build_sketch_dicts(shrinkage_factor):
    # TODO wire offset is wrong for anode or gate
    # Cathode is also shifted
    ptfe_sketches = {}
    ptfe_sketches = TopScreenInsulatingFrame(ptfe_sketches)
    ptfe_sketches = AnodeInsulatingFrame(ptfe_sketches)
    ptfe_sketches = CopperRingInsulatingFrame(ptfe_sketches)
    ptfe_sketches = PTFEWall(ptfe_sketches)
    ptfe_sketches = BottomStackInsulatingFrame(ptfe_sketches)
    ptfe_sketches = CathodeInsulatingFrame(ptfe_sketches)
    ptfe_sketches = GateInsulatingFrame(ptfe_sketches)
    # FIXME What is this component
    #ptfe_sketches = BottomScreenInsulatingFrame(ptfe_sketches)
    ptfe_sketches = PMTReflectors(ptfe_sketches)
    ptfe_sketches = UnderCathodePTFEPin(ptfe_sketches)

    electrode_sketches = {}
    # So far only needed for field shaping guards
    electrode_sketches['shrinkage_factor'] = shrinkage_factor
    electrode_sketches = AllPMTs(electrode_sketches)
    electrode_sketches = TopScreeningElectrode(electrode_sketches)
    electrode_sketches = AnodeElectrode(electrode_sketches)
    electrode_sketches = GateElectrode(electrode_sketches)
    electrode_sketches = CopperRing(electrode_sketches)
    electrode_sketches = FieldShapingRings(electrode_sketches)
    electrode_sketches = FieldShapingGuard(electrode_sketches)
    electrode_sketches = CathodeElectrode(electrode_sketches)
    electrode_sketches = BottomScreeningElectrode(electrode_sketches)
    electrode_sketches = Bell(electrode_sketches)

    xenon_sketches = {}
    xenon_sketches = LXe(xenon_sketches)
    xenon_sketches = GXe(xenon_sketches)

    manual_mapping = {
        1-1: "GXe",
        2-1: "Bell",
        154-1: "GateInsulatingFrame",
        158-1: "GXe",
        402-1: "GateInsulatingFrame",
        479-1: "LXe",
        548-1: "Bell",
    }

    return ptfe_sketches, electrode_sketches, xenon_sketches, manual_mapping