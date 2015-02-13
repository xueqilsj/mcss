from mcss import Potts

import wx
from numpy import zeros, int8, mean, std, exp, log, arange, array, diff
from re import split
from os import getcwd
from time import strftime, time
import threading

import matplotlib
#matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar
import pylab

class CMulticaPotts(Potts):
    def __init__(self):
        Potts.__init__(self)#must use 'self'
        self.RefreshSamp = 100
        self.Stride = 1
        self.Beta0 = 1.42
        self.ReweightSamp = 20
        self.hPrecision = 1e-48
    def Reset(self, dim0, dim1, Q):
        self.InitParameters(Q)
        self.conf_t = zeros((dim0, dim1), int8)
        self.conf = self.PySetConf(self.conf_t)
        self.E = self.PyGetEnergyRange()#numpy array E(n) 
        self.H = zeros((len(self.E), 1), int) #for self.potts.PyTraj2Hist(ydata, self.potts.H)
        self.TR = []#for trajectory
    
    def RunMulX(self):
        Stride = self.Stride * self.nN
        reweightsamples = self.ReweightSamp * self.nN * (self.nHighEBound - self.nLowEBound + 1)
        for CurSample in range(self.RefreshSamp):
            if self.dh > self.hPrecision:
                self.MetroplisTrialSeBound(Stride)
                self.TR.append(self.nCurHami)
                if (self.uCurMove % reweightsamples) == 0:
                    if self.PyMulReweight(self.TR):
                        self.dh /= 2;
                    self.MulReSet()
                    self.TR = []
            else:
                print "Multicanonical method is finished", self.dh 
                break;
                 
                 
          
class FuncPanel(wx.Panel):
    """ A static panel for Potts
    """
    def __init__(self, parent, ID, mainFrame):
        self.potts = CMulticaPotts()
        wx.Panel.__init__(self, parent, ID)
        self.mainFrame = mainFrame;
        self.label_1 = wx.StaticText(self, -1, "L1,L2,Q:  ", style=wx.ALIGN_CENTRE)
        self.text_ctrl_l1l2q = wx.TextCtrl(self, -1, "12,12,10")
        self.label_stride = wx.StaticText(self, -1, "Stride (*N):")
        self.text_ctrl_stride = wx.TextCtrl(self, -1, "3")
        self.label_refreshSamp = wx.StaticText(self, -1, "Refresh Samples:")
        self.text_ctrl_refreshSamp = wx.TextCtrl(self, -1, "1000")
        self.label_reweightSamp = wx.StaticText(self, -1, "Rew. Samples(*Bound):")
        self.text_ctrl_reweightSamp = wx.TextCtrl(self, -1, "20")
        
        self.label_preci = wx.StaticText(self, -1, "Precision of h:")
        self.text_ctrl_preci = wx.TextCtrl(self, -1, "1.0e-48")
        self.label_flatdegree = wx.StaticText(self, -1, "Flat nonzero:")
        self.text_ctrl_flatdegree = wx.TextCtrl(self, -1, "5")
        self.label_beta0 = wx.StaticText(self, -1, "Beta0:")
        self.text_ctrl_beta0 = wx.TextCtrl(self, -1, "1.42")
        self.label_bound = wx.StaticText(self, -1, "Bound (*N):")
        self.text_ctrl_bound = wx.TextCtrl(self, -1, "-1.82:-0.73")
        
        self.label_range = wx.StaticText(self, -1, "From:To ", style=wx.ALIGN_CENTRE)
        self.text_ctrl_range = wx.TextCtrl(self, -1, "0:")
        
        self.label_ranTemp = wx.StaticText(self, -1, "LowT,highT,dT: ", style=wx.ALIGN_CENTRE)
        self.text_ctrl_ranTemp = wx.TextCtrl(self, -1, "0.6, 0.8, 5e-3")
        
        #self.radio_box_initial = wx.RadioBox(self, -1, "Initial state", choices=["Ground", "Random", "Anti-ground", "Backup file"], majorDimension=4, style=wx.RA_SPECIFY_ROWS)
        self.radio_box_leftplot = wx.RadioBox(self, -1, "Left subplot", choices=["Energy Trajectory", "Energy Histogram", "log g(E)", "Diff log g(E)"], majorDimension=4, style=wx.RA_SPECIFY_ROWS)
        self.button_reset = wx.Button(self, -1, "Reset")
        self.button_start = wx.Button(self, -1, "Start")
        self.button_analysis = wx.Button(self, -1, "Analysis")
        self.__do_layout()
        
        self.build_frame(parent)
        
        self.reset = False;
        self.ploting = False
        self.keeping = False
        self.job = None
        self.record_figure_checked = False   
        self.redraw_timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.on_redraw_timer, self.redraw_timer)        
        self.redraw_timer.Start(400)
        #self.on_initial(None)
        self.on_leftplot(None)   

        # end wxGlade
    def __do_layout(self):
        # begin wxGlade: MyFrame.__do_layout
        grp1_box = wx.StaticBox(self, -1, "Potts model")
        grp1_sizer = wx.StaticBoxSizer(grp1_box, wx.HORIZONTAL)
        grp1_box.SetMinSize((220, 142))
        sizer_model1 = wx.BoxSizer(wx.VERTICAL)
        sizer_model1.Add(self.label_1, 0, wx.ALL, 6)
        sizer_model1.Add(self.label_stride, 0, wx.ALL, 6)
        sizer_model1.Add(self.label_refreshSamp, 0, wx.ALL, 6)
        sizer_model1.Add(self.label_reweightSamp, 0, wx.ALL, 6)
        
        sizer_model2 = wx.BoxSizer(wx.VERTICAL)
        sizer_model2.Add(self.text_ctrl_l1l2q, 0, wx.ALL, 2)
        sizer_model2.Add(self.text_ctrl_stride, 0, wx.ALL, 2)
        sizer_model2.Add(self.text_ctrl_refreshSamp, 0, wx.ALL, 2)
        sizer_model2.Add(self.text_ctrl_reweightSamp, 0, wx.ALL, 2)

        self.text_ctrl_l1l2q.SetMinSize((80, 25))
        self.text_ctrl_stride.SetMinSize((80, 25))
        self.text_ctrl_refreshSamp.SetMinSize((80, 25))
        self.text_ctrl_reweightSamp.SetMinSize((80, 25))
                
        self.text_ctrl_stride.Bind(wx.EVT_TEXT, self.on_stride)
        self.text_ctrl_refreshSamp.Bind(wx.EVT_TEXT, self.on_refreshSamp)
        self.text_ctrl_reweightSamp.Bind(wx.EVT_TEXT, self.on_reweightSamp)
               
        grp1_sizer.Add(sizer_model1, 0, wx.ALL, 2)
        grp1_sizer.Add(sizer_model2, 0, wx.ALL, 2)
           
        """        
        grp2_sizer = wx.BoxSizer(wx.VERTICAL)
        grp2_sizer.Add(self.radio_box_initial, 0, wx.ALL, 1)
        self.radio_box_initial.Bind(wx.EVT_RADIOBOX, self.on_initial)
        self.radio_box_initial.SetSelection(1)
        self.radio_box_initial.SetMinSize((120, 150))
        """

        grpflat_box = wx.StaticBox(self, -1, "Reweighting histogram")
        grpflat_sizer = wx.StaticBoxSizer(grpflat_box, wx.HORIZONTAL)
                
        sizer_flathist1 = wx.BoxSizer(wx.VERTICAL)
        sizer_flathist1.Add(self.label_preci, 0, wx.ALL, 6)
        sizer_flathist1.Add(self.label_flatdegree, 0, wx.ALL, 6)
        sizer_flathist1.Add(self.label_beta0, 0, wx.ALL, 6)
        sizer_flathist1.Add(self.label_bound, 0, wx.ALL, 6)
        sizer_flathist2 = wx.BoxSizer(wx.VERTICAL)
        sizer_flathist2.Add(self.text_ctrl_preci, 0, wx.ALL, 2)
        sizer_flathist2.Add(self.text_ctrl_flatdegree, 0, wx.ALL, 2)
        sizer_flathist2.Add(self.text_ctrl_beta0, 0, wx.ALL, 2)   
        sizer_flathist2.Add(self.text_ctrl_bound, 0, wx.ALL, 2)   
        self.text_ctrl_preci.SetMinSize((80, 25))
        self.text_ctrl_flatdegree.SetMinSize((80, 25))
        self.text_ctrl_beta0.SetMinSize((80, 25))
        self.text_ctrl_bound.SetMinSize((80, 25))
        self.text_ctrl_preci.Bind(wx.EVT_TEXT, self.on_preci)
        self.text_ctrl_flatdegree.Bind(wx.EVT_TEXT, self.on_flatdegree)
        self.text_ctrl_beta0.Bind(wx.EVT_TEXT, self.on_beta0)
        self.text_ctrl_bound.Bind(wx.EVT_TEXT, self.on_bound)
                                
        grpflat_sizer.Add(sizer_flathist1, 0, wx.ALL, 2)
        grpflat_sizer.Add(sizer_flathist2, 0, wx.ALL, 2)      
        
        grpsubplot_sizer = wx.BoxSizer(wx.VERTICAL)
        grpsubplot_sizer.Add(self.radio_box_leftplot, 0, wx.ALL, 2)
        self.radio_box_leftplot.SetMinSize((150, 140))  
        self.radio_box_leftplot.Bind(wx.EVT_RADIOBOX, self.on_leftplot)
        self.radio_box_leftplot.SetSelection(1)
     
             
        grpTraj_box = wx.StaticBox(self, -1, "Trajectory Range")
        grpTraj_sizer = wx.StaticBoxSizer(grpTraj_box, wx.VERTICAL)
        grpTraj_box.SetMinSize((118, 38))      
        grpTraj_sizer.Add(self.label_range, 0, wx.ALL, 1)
        grpTraj_sizer.Add(self.text_ctrl_range, 0, wx.ALL, 1)
        self.text_ctrl_range.SetMinSize((100, 26))
        self.text_ctrl_range.Bind(wx.EVT_TEXT, self.on_range)
        
        grp5_box = wx.StaticBox(self, -1, "Analysis")
        grp5_sizer = wx.StaticBoxSizer(grp5_box, wx.VERTICAL)
        grp5_box.SetMinSize((118, 60))
        grp5_sizer.Add(self.label_ranTemp, 0, wx.ALL, 1)
        grp5_sizer.Add(self.text_ctrl_ranTemp, 0, wx.ALL, 1)
        self.text_ctrl_ranTemp.Bind(wx.EVT_TEXT, self.on_ranTemp)
        self.text_ctrl_ranTemp.SetMinSize((100, 26))
        
        grpTJ_sizer = wx.BoxSizer(wx.VERTICAL)
        grpTJ_sizer.Add(grpTraj_sizer, 0, wx.ALL, 1)
        grpTJ_sizer.Add(grp5_sizer, 0, wx.ALL, 1)
 
        grp6_box = wx.StaticBox(self, -1, "")
        grp6_sizer = wx.StaticBoxSizer(grp6_box, wx.VERTICAL)
        grp6_box.SetMinSize((100, 142))  
        grp6_sizer.Add(self.button_reset, border=4, flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        grp6_sizer.Add(self.button_start, border=4, flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        grp6_sizer.Add(self.button_analysis, border=4, flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        
           
        self.button_reset.Bind(wx.EVT_BUTTON, self.on_reset)
        self.button_start.Bind(wx.EVT_BUTTON, self.on_start)
        self.button_analysis.Bind(wx.EVT_BUTTON, self.on_analysis)
        self.button_start.Disable()  
        self.button_analysis.Disable()  
        self.sizer_body = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer_body.Add(grp1_sizer, border=5, flag=wx.ALL)
        #self.sizer_body.Add(grp2_sizer, border=5, flag=wx.ALL)
        self.sizer_body.Add(grpflat_sizer, border=5, flag=wx.ALL)
        #self.sizer_body.AddSpacer(24)
        self.sizer_body.Add(grpsubplot_sizer, border=5, flag=wx.ALL)
        self.sizer_body.Add(grpTJ_sizer, border=5, flag=wx.ALL)
        self.sizer_body.Add(grp6_sizer, border=5, flag=wx.ALL)
        self.SetSizer(self.sizer_body)
        self.sizer_body.Fit(self)
    def build_frame(self, mainPanel):
        self.fig = Figure((8.0, 3.5), dpi=100)
        self.axes_left = self.fig.add_subplot(121)
        self.axes_right = self.fig.add_subplot(122)
        self.axes_left.set_axis_bgcolor('blue')
        self.axes_right.set_axis_bgcolor('blue')
        self.fig.subplots_adjust(bottom=0.14, left=0.11, right=0.96)
        #self.fig.subplots_adjust(top=0.85)
        self.canvas = FigCanvas(mainPanel, -1, self.fig)

        # Create the navigation toolbar, tied to the canvas
        self.toolbar = NavigationToolbar(self.canvas)
        
        # Layout with box sizers
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        #self.vbox.AddSpacer(10)        
        self.vbox.Add(self, 0, flag=wx.ALL | wx.EXPAND)
        mainPanel.SetSizer(self.vbox)
        self.vbox.Fit(self.mainFrame)
    """Event functions"""
    def on_reset(self, event):
        self.reset = False;
        self.button_start.Disable()
        l1l2q = split('[,;]?', self.text_ctrl_l1l2q.GetValue())
        #print l1l2q, bau, self.range 
        if len(l1l2q) != 3:
            self.mainFrame.refresh_status_message("Unknown input:" + self.text_ctrl_l1l2q.GetValue())
            return False
        
        self.potts.Reset(int(l1l2q[0]), int(l1l2q[1]), int(l1l2q[2]))
        self.on_range(None)
        self.on_refreshSamp(None)
        self.on_reweightSamp(None)
        self.on_stride(None)
        self.on_preci(None)
        self.on_flatdegree(None)
        self.on_beta0(None)
        self.on_ranTemp(None)
        self.on_bound(None)     
        self.potts.Beta, self.potts.Se, self.potts.gn, self.potts.hist = self.potts.PySetMul(0.9)
        for i in arange(len(self.potts.Beta)):
            self.potts.Beta[i] = self.potts.Beta0 ;
        self.potts.MulReSet()

        self.potts.WarmupToBound();        
        self.potts.Gauge() 
        self.button_start.Enable()
                
        self.axes_left.clear()
        self.axes_left.text(0.1, 0.85, "Q-state Potts model", color="blue", size=15, transform=self.axes_left.transAxes)
        self.axes_left.text(0.1, 0.6, "$H=-\sum_{<i,j>}\delta(q_i,q_j)$", color="red", size=14, transform=self.axes_left.transAxes)
        self.axes_left.text(0.1, 0.4, "$Z=\sum_Eg(E)e^{-S(E)}$", color="red", size=14, transform=self.axes_left.transAxes)
        self.axes_left.text(0.1, 0.2, "$\\bar Hist^n(E)=max[h,Hist^n(E)]$", color="red", size=14, transform=self.axes_left.transAxes)
        
        self.axes_left.set_axis_bgcolor('black')

        # for configuration plotting
        self.axes_right.clear()
        self.axes_right.set_title(r'Potts Model $N=%d\times%d$,$Q=%d$' % (self.potts.nDim0, self.potts.nDim1, ord(self.potts.cQ)), size=10)
        self.axes_right.set_xlabel('L2', fontsize=9)
        self.axes_right.set_ylabel('L1', fontsize=9)
        self.axes_right.set_axis_bgcolor('black')
        self.axes_right.grid(False)
        self.axes_right.set_aspect("auto", adjustable="box")
        self.im = self.axes_right.imshow(self.potts.conf, interpolation='nearest', origin='lower', animated=True)
        pylab.setp(self.axes_right.get_xticklabels(), fontsize=9)
        pylab.setp(self.axes_right.get_yticklabels(), fontsize=9)
        self.rightplotdirty = True;
        """for tick in self.axes_right.yaxis.get_major_ticks():
            tick.label1On = False
            tick.label2On = True
            #tick.label2.set_color('green')
        """
        self.canvas.draw()  
        
        self.reset = True;
        return True
                 
    def on_start(self, event):
        if not self.keeping:
            self.keeping = True;
            label = "Pause"
            self.button_reset.Disable()
            self.button_analysis.Disable()
            self.text_ctrl_l1l2q.Disable()
            self.text_ctrl_bound.Disable()
            self.starting_time = time()#clock() takes negative when it lasts long  
        else:
            self.keeping = False;
            label = "Start"
            self.button_reset.Enable()
            self.button_analysis.Enable()
            self.text_ctrl_l1l2q.Enable()
            self.text_ctrl_bound.Enable()

        self.button_start.SetLabel(label)
        
    
    def on_redraw_timer(self, event):
        if self.ploting:
            return
        if not self.keeping:
            return
        
        if  self.job == None:
            self.job = threading.Thread(target=self.potts.RunMulX)
            self.job.start()
        else:
            if not self.job.isAlive():
                self.ploting = True
                self.draw_figure()
                self.mainFrame.refresh_status_message("%lg seconds" % (time() - self.starting_time))
                self.job = threading.Thread(target=self.potts.RunMulX)
                self.job.start()
                self.ploting = False
                           
    def draw_figure(self):
        """ Redraws the figure
        """
            
        self.axes_left.clear()
        self.axes_left.set_title('Multicanonical method $h_0=%g$' % self.potts.dh, size=10)
        self.axes_left.grid(True, color='gray')
        self.axes_left.set_axis_bgcolor('black')
        #self.axes_left.lines.remove(self.axes_left_line)
        if self.LeftPlot == 'ETraj':
            if not self.get_data(self.potts.TR):
                return         
            self.axes_left.set_xlabel('Samples', fontsize=9)
            self.axes_left.set_ylabel('Energy', fontsize=9)
            self.axes_left.plot(self.tdata, self.ydata)
            for label in self.axes_left.xaxis.get_ticklabels():
            # label is a Text instance
                label.set_rotation(20)
        elif self.LeftPlot == 'EHist':
            if not self.get_data(self.potts.TR):
                return
            self.axes_left.set_xlabel('Energy', fontsize=9)
            self.axes_left.set_ylabel('Histogram', fontsize=9)
            e0, e1, Erange = self.potts.PyETraj2Hist(self.ydata , self.potts.H)
            #self.axes_left.plot(arange(e0, e1 + 1) / float(self.potts.nN), Erange[e0 + self.potts.nEn - 1:e1 + self.potts.nEn], '.', markersize=3)
            self.axes_left.plot(arange(self.potts.nLowEBound, self.potts. nHighEBound + 1) / float(self.potts.nN), Erange[(self.potts.nLowEBound + self.potts.nEn - 1):(self.potts.nHighEBound + self.potts.nEn)] , 'g.', markersize=2)
        elif self.LeftPlot == 'LogGE':
            self.axes_left.set_xlabel('Energy/N', fontsize=9)
            self.axes_left.set_ylabel(r'$log$ $g(E)$', fontsize=9)
            #self.axes_left.plot(arange(-self.potts.nEn + 1, 1), self.potts.LnGe , 'r.', markersize=2)
            self.axes_left.plot(arange(self.potts.nLowEBound, self.potts. nHighEBound + 1) / float(self.potts.nN), self.potts.Se[(self.potts.nLowEBound + self.potts.nEn - 1):(self.potts.nHighEBound + self.potts.nEn)] , 'r.', markersize=2)
        elif self.LeftPlot == 'DLogGE':
            self.axes_left.set_xlabel('Energy/N', fontsize=9)
            self.axes_left.set_ylabel(r'$Dlog$ $g(E)$', fontsize=9)
            #self.axes_left.plot(arange(-self.potts.nEn + 1, 1), self.potts.LnGe , 'r.', markersize=2)
            self.axes_left.plot(arange(self.potts.nLowEBound, self.potts.nHighEBound) / float(self.potts.nN), self.potts.Beta[(self.potts.nLowEBound + self.potts.nEn - 1):(self.potts.nHighEBound + self.potts.nEn - 1)] , 'c.', markersize=2)

 


        pylab.setp(self.axes_left.get_xticklabels(), fontsize=9)
        pylab.setp(self.axes_left.get_yticklabels(), fontsize=9)
        #plotting configuration  
        if self.rightplotdirty:
            self.axes_right.clear()
            self.axes_right.set_title(r'Potts Model $N=%d\times%d$,$Q=%d$' % (self.potts.nDim0, self.potts.nDim1, ord(self.potts.cQ)), size=10)
            self.axes_right.set_xlabel('L2', fontsize=9)
            self.axes_right.set_ylabel('L1', fontsize=9)
            self.axes_right.set_axis_bgcolor('black')
            self.axes_right.grid(False)
            self.axes_right.set_aspect("auto", adjustable="box")
            self.im = self.axes_right.imshow(self.potts.conf, interpolation='nearest', origin='lower', animated=True)
            pylab.setp(self.axes_right.get_xticklabels(), fontsize=9)
            pylab.setp(self.axes_right.get_yticklabels(), fontsize=9)
            self.rightplotdirty = False
        else: 
            #if self.potts.nInitconf == -1:      
            #    self.axes_right.imshow(self.potts.conf, interpolation='nearest', origin='lower', animated=True)
            self.im.set_array(self.potts.conf) #not suitable for ground state
 
        self.canvas.draw()  
        self.record_figure()
    def get_data(self, TR):
        if len(self.range) == 1:
            return False
        try:
            f = int(self.range[0])
        except:
            self.mainFrame.refresh_status_message("Unknown range:" + ":".join(self.range))
            return False
        
        if self.range[1] == '':
            self.ydata = TR[f:]
            self.tdata = range(len(TR) - len(self.ydata), len(TR));
        else:
            try:
                to = int(self.range[1])
            except:
                self.mainFrame.refresh_status_message("Unknown range:" + ":".join(self.range))
                return False

            self.ydata = TR[f:to]
            x0 = len(TR) - len(self.ydata) if f < 0 else f
            self.tdata = range(x0, x0 + len(self.ydata))
        if len(self.ydata) == 0:
            return False
        return True    
    def exact_range(self, Arr):
        eArr = array(Arr[4:], copy=True)#trim 4 unaccessible energy levels
        eArr[0] = Arr[0]; eArr[1] = Arr[4]
        return eArr
    def on_analysis(self, event):
        #if self.potts.Lnf > self.potts.LnfPrecision:
        #return
        if len(self.ranTemp) != 3:
            self.mainFrame.refresh_status_message("Unknown range of temperature:" + ",".join(self.ranTemp))
            return
        try:
            LowT = float(self.ranTemp[0])
            HighT = float(self.ranTemp[1])
            dt = float(self.ranTemp[2])
        except:
            self.mainFrame.refresh_status_message("Unknown range of temperature:" + ",".join(self.ranTemp))
         
            
        #eE = self.exact_range(self.potts.E)
        #eGeN = self.exact_range(self.potts.GeN)
        #eLnGe = self.exact_range(self.potts.LnGe)
        eE = self.potts.E
        eLnGe = self.potts.Se

        tlist = arange(LowT, HighT, dt).tolist()
        UList = [];FList = []; SList = [];CList = []
        
        lowbi = self.potts.nLowEBound + self.potts.nEn - 1
        highbi = self.potts.nHighEBound + self.potts.nEn - 1
        Level = highbi - lowbi
        for t in tlist:
            e2_l = 0; e_l = 0; z_l = 0
            l = eLnGe[Level / 2] - eE[Level / 2] / t
            for n in range (lowbi, highbi + 1):
                x = eLnGe[n] - eE[n] / t - l
                try:
                    p = exp(x)
                except OverflowError:
                    print "error: exp(x),x=", x, "l=", l
                    return
                z_l += p
                e_l += eE[n] * p
                e2_l += eE[n] * eE[n] * p;
            U = e_l / z_l
            F = -t * (log(z_l) + l)
            C = (e2_l / z_l - U * U) / (t * t)
            S = (U - F) / t
            UList.append(U / self.potts.nN)
            FList.append(F / self.potts.nN)
            CList.append(C / self.potts.nN)
            SList.append(S / self.potts.nN)
        self.axes_left.clear()
        self.axes_left.plot(tlist, UList, 'g.-', label="U")
        self.axes_left.set_xlabel("T")
        self.axes_left.set_title("U/N vs. T, $log(f)=%g$" % self.potts.Lnf, size=10)
        self.axes_left.grid(True, color='gray')
        pylab.setp(self.axes_left.get_xticklabels(), fontsize=9)
        pylab.setp(self.axes_left.get_yticklabels(), fontsize=9)
        
        self.axes_right.clear()
        self.axes_right.set_aspect("auto", adjustable="box", anchor='SW')
        self.axes_right.plot(tlist, CList, 'g.-', label="F")
        self.axes_right.set_xlabel("T")
        self.axes_right.set_title("C/N vs. T $log(f)=%g$" % self.potts.Lnf, size=10)
        self.axes_right.grid(True, color='gray')
        pylab.setp(self.axes_right.get_xticklabels(), fontsize=9)
        pylab.setp(self.axes_right.get_yticklabels(), fontsize=9)
        
        self.canvas.draw()  
        self.record_figure() 
        self.rightplotdirty = True;
        
    def on_ranTemp(self, event):
        self.ranTemp = split('[,;:]?', self.text_ctrl_ranTemp.GetValue())
    def on_initial(self, event):
        self.curchoice = self.radio_box_initial.GetSelection()
        if self.curchoice == 0:
            self.initConf = "-1"
        elif self.curchoice == 1:
            self.initConf = "0"
        elif self.curchoice == 2:
            self.initConf = "1"
        else:
            file_choices = "Conf (*.txt)|*.txt"
            self.ConfPath = ""
            dlg = wx.FileDialog(
                self,
                message="Open configurational file...",
                defaultDir=getcwd(),
                defaultFile="",
                wildcard=file_choices,
                style=wx.OPEN)
        
            if dlg.ShowModal() == wx.ID_OK:
                self.ConfPath = dlg.GetPath()
                self.initConf = "2"
        #event.Skip()
  
    def on_leftplot(self, event):
        curchoice = self.radio_box_leftplot.GetSelection()
        if curchoice == 0:
            self.LeftPlot = "ETraj"
        elif curchoice == 1:
            self.LeftPlot = "EHist"
        elif curchoice == 2:
            self.LeftPlot = "LogGE"
        elif curchoice == 3:
            self.LeftPlot = "DLogGE"
        #event.Skip()
    def on_preci(self, event):
        try:
            f = float(self.text_ctrl_preci.GetValue())
        except:
            self.mainFrame.refresh_status_message("Unknown Precision:", self.text_ctrl_preci.GetValue())
            return
        if f < 0:
            self.mainFrame.refresh_status_message("Too small value %g" % f + " and reset to 1e-8")
            f = 1e-8
            return
        self.potts.hPrecision = f;
        return
    def on_flatdegree(self, event):
        try:
            f = int(self.text_ctrl_flatdegree.GetValue())
        except:
            self.mainFrame.refresh_status_message("Unknown Flat degree:", self.text_ctrl_flatdegree.GetValue())
            return
        if f < 1:
            self.mainFrame.refresh_status_message("Illegal value %g" % f + " and reset to 1")
            f = 1
            return
        self.potts.nFlatNonzero = f;   
        return
    def on_beta0(self, event):
        try:
            f = float(self.text_ctrl_beta0.GetValue())
        except:
            self.mainFrame.refresh_status_message("Unknown Flat beta0:", self.text_ctrl_beta0.GetValue())
            return
        if f < 0:
            self.mainFrame.refresh_status_message("Illegal value %g" % f + " and reset to 1.42")
            f = 1.42
            return
        self.potts.Beta0 = f;  
        return
    def on_bound(self, event):
        bound = split('[,;:]?', self.text_ctrl_bound.GetValue())
        if len(bound) == 1:
            self.mainFrame.refresh_status_message("Unknown bound:", self.text_ctrl_bound.GetValue())
            return
        try :
            lb = float(bound[0])
        except:
            self.mainFrame.refresh_status_message("Unknown bound:", self.text_ctrl_bound.GetValue())
            return

        lb6 = int(lb * self.potts.nN)
        if lb > 0 or lb6 < -2 * self.potts.nN - 5: # lb >= 6
            self.mainFrame.refresh_status_message("Illegal lower bound %g" % lb)                
            return
        
        if bound[1] == "":
            ub = 2;
        else: 
            try :
                ub = float(bound[1])
            except:
                self.mainFrame.refresh_status_message("Unknown bound:", self.text_ctrl_bound.GetValue())
                return
            if ub > 0 or ub < -2:
                self.mainFrame.refresh_status_message("Illegal upper bound %g" % ub)
                return
        nLowEBound = int(lb * self.potts.nN)
        nHighEBound = int(ub * self.potts.nN)
        if nLowEBound >= nHighEBound:
            self.mainFrame.refresh_status_message("nLowEBound %d great than upper nHighEBound %d" % (nLowEBound, nHighEBound))
            return
        self.potts.nLowEBound = nLowEBound
        self.potts.nHighEBound = nHighEBound
        print "Bound:", nLowEBound, nHighEBound
    def on_range(self, event):
        self.range = split('[,;:]?', self.text_ctrl_range.GetValue())
    def on_refreshSamp(self, event):
        try:
            f = int(self.text_ctrl_refreshSamp.GetValue())
        except:
            self.mainFrame.refresh_status_message("Unknown Refresh Samples:", self.text_ctrl_refreshSamp.GetValue())
            return
        if f < 1:
            self.mainFrame.refresh_status_message("Too small value %d" % f + " and reset to 1")
            f = 1
            return
        self.potts.RefreshSamp = f;    
    def on_reweightSamp(self, event):
        try:
            f = int(self.text_ctrl_reweightSamp.GetValue())
        except:
            self.mainFrame.refresh_status_message("Unknown Reweight Samples:", self.text_ctrl_reweightSamp.GetValue())
            return
        if f < 1:
            self.mainFrame.refresh_status_message("Too small value %d" % f + " and reset to 1")
            f = 1
            return
        self.potts.ReweightSamp = f;   
    def on_stride(self, event):
        try:
            f = int(self.text_ctrl_stride.GetValue())
        except:
            self.mainFrame.refresh_status_message("Unknown Stride:", self.text_ctrl_stride.GetValue())
            return
        if f < 1:
            self.mainFrame.refresh_status_message("Too small value %d" % f + " and reset to 1")
            f = 1
            return
        self.potts.Stride = f;    
       
    """for Mainframe calling functions"""
    def message_dlg(self, titile, msg):
            dlg = wx.MessageDialog(self, msg, titile, wx.OK | wx.ICON_INFORMATION)
            dlg.ShowModal() 
            dlg.Destroy()
            
    def save_traj(self):
        if self.keeping: 
            self.message_dlg("Error", "Stop runing first!")
            return
        if not self.reset: 
            self.message_dlg("Error", "No resetting!")
            return
        
        file_choices = "lnGe_GeN (*.txt)|*.txt"
        trajfile = "wlxp%d_%d_%d~%lg_%lg_%d_%d_%dLnGe_GeN.txt" % (self.potts.nDim0, self.potts.nDim1, ord(self.potts.cQ), \
            self.potts.LnfPrecision, self.potts.FlatRate, self.potts.Stride, \
            self.potts.nLowEBound, self.potts.nHighEBound);
                        
        dlg = wx.FileDialog(self, message="Save trajectory file...",
                defaultDir=getcwd(), defaultFile=trajfile, wildcard=file_choices,
                style=wx.SAVE)
                
        if dlg.ShowModal() == wx.ID_OK:
                trajfile = dlg.GetPath()
                f = open(trajfile, 'w')
                if f:
                    for n in range(len(self.potts.LnGe)):
                        f.write(str(self.potts.LnGe[n]) + "\t" + str(self.potts.GeN[n]) + "\n")

    def save_conf(self):
        if self.keeping: 
            self.message_dlg("Error", "Stop runing first!")
            return
        if not self.reset: 
            self.message_dlg("Error", "No resetting!")
            return 
        
        file_choices = "Conf (*.conf)|*.conf"
        if self.using_wolff:
            trajfile = "gcewfp%d_%d_%d~%d~%.15lg_%.15lg_%.15lg_%ld_%ld.conf" % (self.potts.nDim0, self.potts.nDim1, ord(self.potts.cQ), self.potts.nInitconf, \
            self.potts.dBeta, self.potts.dAlpha, self.potts.dUn, \
            self.potts.Stride, self.potts.RefreshSamp);
        else:
            trajfile = "gce%d_%d_%d~%d~%.15lg_%.15lg_%.15lg_%ld_%ld.conf" % (self.potts.nDim0, self.potts.nDim1, ord(self.potts.cQ), self.potts.nInitconf, \
            self.potts.dBeta, self.potts.dAlpha, self.potts.dUn, \
            self.potts.Stride, self.potts.RefreshSamp);
                        
        dlg = wx.FileDialog(self, message="Save configuration file...",
                defaultDir=getcwd(), defaultFile=trajfile, wildcard=file_choices,
                style=wx.SAVE)
                
        if dlg.ShowModal() == wx.ID_OK:
                trajfile = dlg.GetPath()
                print trajfile
                self.potts.DumpConf(trajfile)
    def record_figure(self):
        if self.record_figure_checked:
            figfile = "mulp%d_%d_%d~%d~%.15lg_%.15lg_%ld_%ld_%ld_%ldplot%s.png" % (self.potts.nDim0, self.potts.nDim1, ord(self.potts.cQ), self.potts.nInitconf, \
            self.potts.Beta0, self.potts.hPrecision, self.potts.ReweightSamp, \
            self.potts.Stride, self.potts.RefreshSamp, self.potts.uCurMove, strftime("%Y-%m-%d@%H~%M~%S"));
            self.fig.savefig(figfile) 

            print "save figure to:", figfile
        else:
            pass
