
from pylab import *
import tkinter as tk
import control as ctrl
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure


window = tk.Tk()
f = Figure()

class Spec:
    wc = tk.DoubleVar(window,value=0.2)
    pm = tk.DoubleVar(window,value=60)
    A = tk.DoubleVar(window,value=1)
    B = tk.DoubleVar(window,value=5)
    e = tk.DoubleVar(window,value=0.1)

    def Conf_Spec():
        Spec.wc = double(ewc.get())
        Spec.pm = double(epm.get())
        Spec.A = double(eA.get())
        Spec.B = double(eB.get())
        Spec.e = double(ee.get())

class G:
    tau1 = tk.DoubleVar(window,value=1)
    tau2 = tk.DoubleVar(window,value=5)
    tau3 = tk.DoubleVar(window,value=10)
    DCgain = tk.DoubleVar(window,value=10)

    NUM = [1]
    DEN = [1, 1, 1]
    sysG = ctrl.tf(NUM, DEN)

    def Conf_G( ):
        G.tau1 = float(etau1.get())
        G.tau2 = float(etau2.get())
        G.tau3 = float(etau3.get())
        G.DCgain = float(eDCgain.get())

        G.NUM = [G.DCgain]
        G.DEN = [G.tau1 * G.tau2 * G.tau3, G.tau1 * G.tau2 + G.tau2 * G.tau3 + G.tau1 * G.tau3,
                 G.tau1 + G.tau2 + G.tau3, 1]
        G.sysG = ctrl.tf(G.NUM, G.DEN)
class R:
    Case=tk.IntVar(window,value=0)
    ur = 1
    urOffset=1 # used in case c to adjust the DC gain of R
    R2NUM = [1]
    R2DEN = [1, 1, 1]
    sysR1= ctrl.tf([1],[1])
    sysR2= ctrl.tf(R2NUM, R2DEN)
    sysR = ctrl.tf([1],[1])

    def CaseA(sysg,wctar):
        R.ur = ((Spec.A + Spec.B) / Spec.e - 1) * 1.3 / G.DCgain
        R.sysR1 = ctrl.tf([R.ur], [1])
        tauAlow = ctrl.dcgain(G.sysG*R.sysR1) / (wctar*1.5)
        tauAhigh = sqrt(prod(-ctrl.pole(G.sysG*R.sysR1) ** -1) / tauAlow)
        R.R2NUM=G.DEN
        R.R2DEN=[tauAhigh * tauAhigh * tauAlow, 2 * tauAhigh * tauAlow, 2 * tauAhigh + tauAlow, 1]
        R.sysR2 = ctrl.tf(R.R2NUM,R.R2DEN )
        R.sysR=R.sysR1*R.sysR2
        return R.sysR

    def CaseB(sysg,wctar):
        para=50
        R.ur = ((Spec.A + Spec.B) / Spec.e - 1) * 1.3 / G.DCgain
        R.sysR1 = ctrl.tf([R.ur], [1])
        R.R2NUM=[(-ctrl.pole(G.sysG)[0]**-1)*(-ctrl.pole(G.sysG)[1]**-1), (-ctrl.pole(G.sysG)[0]**-1)+(-ctrl.pole(G.sysG)[1]**-1),1]
        R.R2DEN=[(-para*ctrl.pole(G.sysG)[0])**-2, 2*((-para*ctrl.pole(G.sysG)[0])**-1), 1]
        R.sysR2 = ctrl.tf(R.R2NUM, R.R2DEN)
        R.sysR = R.sysR1 * R.sysR2
        return R.sysR

    def CaseC(sysg,wctar):
        PoleCanN=2
        if PoleCanN == 0:
            R.R2NUM = [1]
            R.R2DEN = [1]

        elif PoleCanN == 1:
            R.R2NUM = [-ctrl.pole(G.sysG)[2] ** -1, 1]
            R.R2DEN = [1]

        elif PoleCanN == 2:
            R.R2NUM = [(-ctrl.pole(G.sysG)[1] ** -1) * (-ctrl.pole(G.sysG)[2] ** -1),
                       (-ctrl.pole(G.sysG)[1] ** -1) + (-ctrl.pole(G.sysG)[2] ** -1), 1]
            R.R2DEN = [-ctrl.pole(G.sysG)[0] ** -1, 1]

        R.sysR1 = ctrl.tf([R.ur], [1, 0])
        R.sysR2 = ctrl.tf(R.R2NUM, R.R2DEN)
        R.sysR = R.sysR1 * R.sysR2

        return R.sysR

    def urCheck(sys,wctar,pmtar,omega):
        thres=0.5
        R.urOffset = 1
        gm, pm, wg, wp = ctrl.margin(sys)
        mag, phaserad, w = ctrl.bode(sys, omega)
        phasedeg = phaserad * 180 / pi
        if pm < pmtar:
            for i in range(len(omega)):
                if phasedeg[i] < pmtar-180+thres and phasedeg[i] > pmtar-180-thres and mag[i]>1:
                    R.urOffset=mag[i]**-1
                    break
        return sys * R.urOffset
    def Conf_R():
        R.Case = int(ecase.get())
        if R.Case == 1:
            R.sysR = R.CaseA(G.sysG, Spec.wc)
        elif R.Case == 2:
            R.sysR = R.CaseB(G.sysG, Spec.wc)
        elif R.Case == 3:
            R.sysR = R.CaseC(G.sysG, Spec.wc)

class L:
    sysL = ctrl.tf([1], [1])
    wc_re = tk.DoubleVar(window)
    pm_re = tk.DoubleVar(window)
    def Bode_Margin(sys, omega):
        w = omega
        mag, phaserad, w = ctrl.bode(sys, w)
        magdb = 20 * log(mag)
        phasedeg = phaserad * 180 / pi
        gm, pm, wg, wp = ctrl.margin(sys)

        a1 = f.add_subplot(221)
        a1.set_xscale("log")
        a1.plot(w, magdb)
        a1.plot(wp, 0, '.y')
        a1.text(wp, 0 + 25, "wc=%.2f" % wp)
        # a1.title = "Bode"
        # a1.ylabel = "Magnitude dB "
        # a1.xlabel = "Angular frequency"

        a2 = f.add_subplot(223)
        a2.set_xscale("log")
        a2.plot(w, phasedeg)
        a2.plot(wp, pm - 180, '.y')
        a2.text(wp, pm - 180 + 25, "phase=%.2f" % (pm - 180))
        # a2.ylabel = "Phase (degree) "
        # a2.xlabel = "Angular frequency"

        a3 = f.add_subplot(122)
        sysF = sys / (1 + sys)
        T = np.arange(0, 30, 0.5)
        T, yout = ctrl.step_response(sysF, T)
        a3.plot(T, yout)
        # a3.title = "step response"
        # a3.xlabel = "time"
        # a3.ylabel = "Magnitude"

        return wp, pm



    def Conf_L( ):
        f.clf()
        Spec.Conf_Spec()
        G.Conf_G()
        R.Conf_R()
        N = 1000
        fr = logspace(-5, 1, N)
        w = 2 * pi * fr
        if R.Case == 3:
            L.sysL = R.urCheck(G.sysG * R.sysR, Spec.wc, Spec.pm, w)
        else :
            L.sysL = G.sysG * R.sysR
        L.wc_re, L.pm_re = L.Bode_Margin(L.sysL, w)

        ewc_re.delete(0, "end")
        ewc_re.insert(0, '%.2f'%L.wc_re)
        epm_re.delete(0, "end")
        epm_re.insert(0,'%.2f'%L.pm_re)
        cv.draw()

# window
window.title('Controller Design')
window.geometry('1080x1600')
# Spec
tk.Label(window, text='Specification',fg='blue', font=('Arial', 12), justify='left', width=21, height=2,)\
    .place(x=15, y=10, anchor='nw')
tk.Label(window, text='critical frequency >', font=('Arial', 12), justify='left', width=21, height=2)\
    .place(x=15, y=40, anchor='nw')
tk.Label(window, text=' marginal phase >', font=('Arial', 12), justify='left', width=21, height=2)\
    .place(x=15, y=70, anchor='nw')
tk.Label(window, text='magnitude of input=', font=('Arial', 12), justify='left', width=21, height=2)\
    .place(x=15, y=100, anchor='nw')
tk.Label(window, text='magnitude of disturbance=', font=('Arial', 12), justify='left', width=21, height=2)\
    .place(x=15, y=130, anchor='nw')
tk.Label(window, text='steady state of error<', font=('Arial', 12), justify='left', width=21, height=2)\
    .place(x=15, y=160, anchor='nw')
ewc = tk.Entry(window, font=('Arial', 10) ,width=10, text=Spec.wc)
ewc.place(x=210, y=50, anchor='nw')
epm = tk.Entry(window, font=('Arial', 10), width=10, text=Spec.pm)
epm.place(x=210, y=80, anchor='nw')
eA = tk.Entry(window, font=('Arial', 10), width=10, text=Spec.A)
eA.place(x=210, y=110, anchor='nw')
eB = tk.Entry(window, font=('Arial', 10),width=10, text=Spec.B)
eB.place(x=210, y=140, anchor='nw')
ee = tk.Entry(window, font=('Arial', 10), width=10, text=Spec.e)
ee.place(x=210, y=170, anchor='nw')
# tk.Button(window, text='Confirm specification', font=('Arial', 12), width=25, height=1, command=Spec.Conf_Spec)\
#     .place(x=50, y=220, anchor='nw')
# G
tk.Label(window, text='Parameters', fg='blue', font=('Arial', 12), justify='left', width=21, height=2)\
    .place(x=315, y=10, anchor='nw')
tk.Label(window, text='t1=', font=('Arial', 12), justify='left', width=21, height=2)\
    .place(x=315, y=40, anchor='nw')
tk.Label(window, text='t2=', font=('Arial', 12), justify='left', width=21, height=2)\
    .place(x=315, y=70, anchor='nw')
tk.Label(window, text='t3=', font=('Arial', 12), justify='left', width=21, height=2)\
    .place(x=315, y=100, anchor='nw')
tk.Label(window, text='DC gain=', font=('Arial', 12), justify='left', width=21, height=2)\
    .place(x=315, y=130, anchor='nw')
etau1 = tk.Entry(window, font=('Arial', 10),width=10, text=G.tau1)
etau1.place(x=510, y=50, anchor='nw')
etau2 = tk.Entry(window, font=('Arial', 10),width=10, text=G.tau2)
etau2.place(x=510, y=80, anchor='nw')
etau3 = tk.Entry(window, font=('Arial', 10),width=10, text=G.tau3)
etau3.place(x=510, y=110, anchor='nw')
eDCgain = tk.Entry(window, font=('Arial', 10),width=10, text=G.DCgain)
eDCgain.place(x=510, y=140, anchor='nw')

# tk.Button(window, text='Confirm object', font=('Arial', 12), width=25, height=1, command=G.Conf_G)\
#     .place(x=350, y=220, anchor='nw')
# R
tk.Label(window, text='Controller design', fg='blue', font=('Arial', 12), justify='left', width=21, height=2,)\
    .place(x=615, y=10, anchor='nw')
tk.Label(window, text='Select a case(1/2/3):', font=('Arial', 12), justify='left', width=21, height=2,)\
    .place(x=615, y=40, anchor='nw')
ecase = tk.Entry(window, font=('Arial', 10),width=10,text=R.Case)
ecase.place(x=810, y=50, anchor='nw')

# tk.Button(window, text='Confirm', font=('Arial', 12), width=25, height=1, command=R.Conf_R)\
#     .place(x=650, y=220, anchor='nw')

# result
tk.Label(window, text='Result', fg='blue', font=('Arial', 12), justify='left', width=21, height=2,)\
    .place(x=15, y=300, anchor='nw')
tk.Label(window, text='critical frequency =', font=('Arial', 12), justify='left', width=21, height=2)\
    .place(x=145, y=300, anchor='nw')
tk.Label(window, text=' marginal phase =', font=('Arial', 12), justify='left', width=21, height=2)\
    .place(x=385, y=300, anchor='nw')
ewc_re = tk.Entry(window, font=('Arial', 10),width=10)
ewc_re.place(x=315, y=310, anchor='nw')
epm_re = tk.Entry(window, font=('Arial', 10),width=10)
epm_re.place(x=555, y=310, anchor='nw')

tk.Button(window, text='Confirm', font=('Arial', 14), width=35, height=1, command=L.Conf_L)\
    .place(x=310, y=260, anchor='nw')

cv = FigureCanvasTkAgg(f,window)
cv.get_tk_widget().place(x=145, y=380, anchor='nw')
cv.draw()

def main():
    window.mainloop()
    print(G.sysG)
    print(R.sysR)
    print(L.sysL)

if __name__ == "__main__":
    main()
