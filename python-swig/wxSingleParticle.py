from mcss import SParticle, PyJackknife_MeanStdErrTau_D

import wx
from numpy import zeros, int8, mean, std, arange, array
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

class CSingleParticle(SParticle):
    def __init__(self):
        SParticle.__init__(self)#must use 'self'
        self.RefreshSamp = 100
        self.Stride = 1
        self.Potential = None;
    def Reset(self, initX, off, T):
        self.InitParameters(initX, off, T)
        self.TR = []#for trajectory
        self.TRX = []#for trajectory
                
    def GetXY(self, step):
        if self.Potential == "DB":
            return self.PyGetXY_Pak(step)
        elif self.Potential == "Mul":
            return self.PyGetXY_MP(step)
        else:
            print "Custom"
    def GetXlim(self):
        if self.Potential == "DB":
            return self.PAK_XL, self.PAK_XH
        elif self.Potential == "Mul":
            return self.MP_XL, self.MP_XH
        else:
            print "Custom"
    def GetYlim(self):
        if self.Potential == "DB":
            return - 4.08304839, 4.73578521
        elif self.Potential == "Mul":
            return 0, 8     
        else:
            print "Custom"
            return None
    def Guage(self):
        if self.Potential == "DB":
            return self.Gauge_Pak()
        elif self.Potential == "Mul":
            return self.Gauge_MP()
        else:
            print "Custom"

    def RunMetroplisTrial(self):
        if self.Potential == "DB":
            #self.Gauge_Pak()#initp = 0, offset = 0.002, temp = 1 / 5.0
            for CurSample in range(self.RefreshSamp):
                self.MetroplisTrial_Pak(self.Stride)
                self.TR.append(self.GetArg(self.AI_HAMI))
                self.TRX.append(self.dX)
        elif self.Potential == "Mul":
            #self.Gauge_MP()
            for CurSample in range(self.RefreshSamp):
                self.MetroplisTrial_MP(self.Stride)
                self.TR.append(self.GetArg(self.AI_HAMI))
                self.TRX.append(self.dX)

          
class FuncPanel(wx.Panel):
    """ A static panel for sparticle
    """
    def __init__(self, parent, ID, mainFrame):
        self.sparticle = CSingleParticle()
        wx.Panel.__init__(self, parent, ID)
        self.mainFrame = mainFrame;
      
        self.label_1 = wx.StaticText(self, -1, "Initial position:", style=wx.ALIGN_CENTRE)
        self.text_ctrl_initial = wx.TextCtrl(self, -1, "0")
        self.label_2 = wx.StaticText(self, -1, "Offset position:")
        self.text_ctrl_offset = wx.TextCtrl(self, -1, "0.1")
        self.label_3 = wx.StaticText(self, -1, "Temperature:")
        self.text_ctrl_temperature = wx.TextCtrl(self, -1, "1.5")
        
        self.label_stride = wx.StaticText(self, -1, "Stride:")
        self.text_ctrl_stride = wx.TextCtrl(self, -1, "1")
        self.label_refreshSamp = wx.StaticText(self, -1, "Refresh Samples:")
        self.text_ctrl_refreshSamp = wx.TextCtrl(self, -1, "100")

        self.label_range = wx.StaticText(self, -1, "From:To ", style=wx.ALIGN_CENTRE)
        self.text_ctrl_range = wx.TextCtrl(self, -1, "-2000:")
    

        self.radio_box_potential = wx.RadioBox(self, -1, "Potential function", choices=["Double wells", "Multiple wells"], majorDimension=3, style=wx.RA_SPECIFY_ROWS)
        self.radio_box_upperplot = wx.RadioBox(self, -1, "Upper subplot", choices=["X positon trajectory", "X positon histogram", "Potential F(x) trajectory", "Potential F(x) histogram"], majorDimension=4, style=wx.RA_SPECIFY_ROWS)
        self.button_start = wx.Button(self, -1, "Start")
        self.button_analysis = wx.Button(self, -1, "Analysis")
        self.__do_layout()
        
        self.build_frame(parent)
        
        self.reset = False;
        self.ploting = False
        self.keeping = False
        self.job = None
        self.record_figure_checked = False   
        self.transX = matplotlib.transforms.blended_transform_factory(self.axes_upper.transData, self.axes_upper.transAxes)    
        self.transY = matplotlib.transforms.blended_transform_factory(self.axes_upper.transAxes, self.axes_upper.transData)

        self.redraw_timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.on_redraw_timer, self.redraw_timer)        
        self.redraw_timer.Start(200)
        self.on_potential(None)

        # end wxGlade
    def __do_layout(self):
        # begin wxGlade: MyFrame.__do_layout
        grp1_sizer = wx.BoxSizer(wx.VERTICAL)
        grp1_sizer.Add(self.radio_box_potential, 0, wx.ALL, 2)
        self.radio_box_potential.Bind(wx.EVT_RADIOBOX, self.on_potential)
        self.radio_box_potential.SetSelection(1)
        self.radio_box_potential.SetMinSize((220, 70))

        
        grp2_box = wx.StaticBox(self, -1, "Sparticle model")
        grp2_box.SetMinSize((220, 100))
        grp2_sizer = wx.StaticBoxSizer(grp2_box, wx.HORIZONTAL)

        sizer_model1 = wx.BoxSizer(wx.VERTICAL)
        sizer_model1.Add(self.label_1, 0, wx.ALL, 6)
        sizer_model1.Add(self.label_2, 0, wx.ALL, 6)
        sizer_model1.Add(self.label_3, 0, wx.ALL, 6)
        sizer_model2 = wx.BoxSizer(wx.VERTICAL)
        sizer_model2.Add(self.text_ctrl_initial, 0, wx.ALL, 1)
        sizer_model2.Add(self.text_ctrl_offset, 0, wx.ALL, 1)
        sizer_model2.Add(self.text_ctrl_temperature, 0, wx.ALL, 1)
        self.text_ctrl_temperature.Bind(wx.EVT_TEXT, self.on_temperature)
        
        grp2_sizer.Add(sizer_model1, 0, wx.ALL, 2)
        grp2_sizer.Add(sizer_model2, 0, wx.ALL, 2)
                   
        grp3_box = wx.StaticBox(self, -1, "Sampling")
        grp3_box.SetMinSize((220, 70))
        grp3_sizer = wx.StaticBoxSizer(grp3_box, wx.HORIZONTAL)
        
        sizer_sample1 = wx.BoxSizer(wx.VERTICAL)
        sizer_sample1.Add(self.label_stride, 0, wx.ALL, 6)
        sizer_sample1.Add(self.label_refreshSamp, 0, wx.ALL, 6)
        sizer_sample2 = wx.BoxSizer(wx.VERTICAL)
        sizer_sample2.Add(self.text_ctrl_stride, 0, wx.ALL, 2)
        sizer_sample2.Add(self.text_ctrl_refreshSamp, 0, wx.ALL, 2)
        self.text_ctrl_stride.Bind(wx.EVT_TEXT, self.on_stride)
        self.text_ctrl_refreshSamp.Bind(wx.EVT_TEXT, self.on_refreshSamp)
        
        grp3_sizer.Add(sizer_sample1, 0, wx.ALL, 2)
        grp3_sizer.Add(sizer_sample2, 0, wx.ALL, 2)

        
        grp4_sizer = wx.BoxSizer(wx.VERTICAL)
        grp4_sizer.Add(self.radio_box_upperplot, 0, wx.ALL, 2)
        self.radio_box_upperplot.SetMinSize((220, 120))
        self.radio_box_upperplot.Bind(wx.EVT_RADIOBOX, self.on_upperplot)
        self.radio_box_upperplot.SetSelection(0)
        #self.on_upperplot(None)        
        grp5_box = wx.StaticBox(self, -1, "Trajectory Range")
        grp5_sizer = wx.StaticBoxSizer(grp5_box, wx.HORIZONTAL)  
        grp5_box.SetMinSize((220, 30))
        grp5_sizer.Add(self.label_range, 0, wx.ALL, 6)
        grp5_sizer.Add(self.text_ctrl_range, 0, wx.ALL, 1)
        self.text_ctrl_range.Bind(wx.EVT_TEXT, self.on_range)
        grp4_sizer.Add(grp5_sizer, 0, wx.ALL, 2)

 
        grp6_box = wx.StaticBox(self, -1, "")
        grp6_sizer = wx.StaticBoxSizer(grp6_box, wx.HORIZONTAL)
        grp6_box.SetMinSize((220, 30))
        grp6_sizer.Add(self.button_analysis, border=4, flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        grp6_sizer.Add(self.button_start, border=4, flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL)
           
        self.button_start.Bind(wx.EVT_BUTTON, self.on_start)
        self.button_analysis.Bind(wx.EVT_BUTTON, self.on_analysis)
        self.button_start.Disable()  
        self.button_analysis.Disable()
        self.sizer_body = wx.BoxSizer(wx.VERTICAL)
        self.sizer_body.Add(grp1_sizer, border=5, flag=wx.ALL)
        self.sizer_body.Add(grp2_sizer, border=5, flag=wx.ALL)
        self.sizer_body.Add(grp3_sizer, border=5, flag=wx.ALL)
        #self.sizer_body.AddSpacer(24)
        self.sizer_body.Add(grp4_sizer, border=5, flag=wx.ALL)
        self.sizer_body.Add(grp6_sizer, border=5, flag=wx.ALL)
        self.SetSizer(self.sizer_body)
        self.sizer_body.Fit(self)
    def build_frame(self, mainPanel):
        self.fig = Figure((4.5, 6.5), dpi=100)
        self.axes_upper = self.fig.add_subplot(211)
        self.axes_lower = self.fig.add_subplot(212)
        self.fig.subplots_adjust(bottom=0.06, top=0.95, left=0.15, right=0.95)
        #self.fig.subplots_adjust(top=0.85)
        self.canvas = FigCanvas(mainPanel, -1, self.fig)
        # Bind the 'pick' event for clicking on one of the bars
        # Create the navigation toolbar, tied to the canvas
        self.toolbar = NavigationToolbar(self.canvas)
        
        # Layout with box sizers
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
                
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        #self.vbox.AddSpacer(10)    
        self.hbox.Add(self.vbox, 0, flag=wx.ALL | wx.EXPAND)
        self.hbox.Add(self, 0, flag=wx.ALL | wx.EXPAND)
        mainPanel.SetSizer(self.hbox)
        self.hbox.Fit(self.mainFrame)

    def on_reset(self, event):
        self.reset = False;
        self.button_start.Disable()
        try:
            initX = float(self.text_ctrl_initial.GetValue())
        except:
            self.mainFrame.refresh_status_message("Unknown initial position:" + self.text_ctrl_initial.GetValue())
            return
        try:
            off = float(self.text_ctrl_offset.GetValue())
        except:
            self.mainFrame.refresh_status_message("Unknown offset position:" + self.text_ctrl_offset.GetValue())
            return
        try:
            T = float(self.text_ctrl_temperature.GetValue())
        except:
            self.mainFrame.refresh_status_message("Unknown inverse temperature:" + self.text_ctrl_temperature.GetValue())
            return
        self.sparticle.Reset(initX, off, T)
        self.on_refreshSamp(None)
        self.on_stride(None)
        self.on_range(None)
        self.on_upperplot(None)

        self.button_start.Enable()
                
        self.axes_upper.clear()
        self.axes_upper.set_axis_bgcolor('black')
        pylab.setp(self.axes_upper.get_xticklabels(), fontsize=9)
        pylab.setp(self.axes_upper.get_yticklabels(), fontsize=9)
        # for configuration plotting
        self.axes_lower.clear()
        self.axes_lower.set_xlabel('X', fontsize=9)
        self.axes_lower.set_ylabel('Potential f(X)', fontsize=9)
        self.axes_lower.set_axis_bgcolor('black')
        self.axes_lower.grid(True, color='gray')
        pylab.setp(self.axes_lower.get_xticklabels(), fontsize=9)
        pylab.setp(self.axes_lower.get_yticklabels(), fontsize=9)
        XY = self.sparticle.GetXY(0.01)
        self.axes_lower.plot(XY[0], XY[1], 'g-');
        xl0, xl1 = self.sparticle.GetXlim()
        self.axes_lower.set_xlim(xl0, xl1)
        
        self.sparticle.Guage()
        self.curXYpoint, = self.axes_lower.plot([self.sparticle.dX], [self.sparticle.GetArg(self.sparticle.AI_HAMI)], 'ro')
        self.canvas.draw()  
        
        self.reset = True;
        return True
                 
    def on_start(self, event):
        if not self.keeping:
            self.keeping = True;
            label = "Pause"
            self.button_analysis.Disable()
            self.radio_box_potential.Disable()
            self.text_ctrl_initial.Disable()    
            self.text_ctrl_offset.Disable()
        
            self.starting_time = time()#clock() takes negative when it lasts long  
        else:
            self.keeping = False;
            label = "Start"
            self.button_analysis.Enable()
            self.radio_box_potential.Enable()
            self.text_ctrl_initial.Enable()    
            self.text_ctrl_offset.Enable()

        self.button_start.SetLabel(label)
        
    
    def on_redraw_timer(self, event):
        if self.ploting:
            return
        if not self.keeping:
            return
        
        if  self.job == None:
            self.job = threading.Thread(target=self.sparticle.RunMetroplisTrial)
            self.job.start()
        else:
            if not self.job.isAlive():
                self.ploting = True
                self.draw_figure()
                self.mainFrame.refresh_status_message("%lg seconds" % (time() - self.starting_time))
                self.job = threading.Thread(target=self.sparticle.RunMetroplisTrial)
                self.job.start()
                self.ploting = False
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
    def draw_figure(self):
        """ Redraws the figure
        """
        self.axes_upper.clear()
        self.axes_upper.set_title(r'Temperature $T=%lg$' % (self.sparticle.GetArg(self.sparticle.AI_TEMPER)), size=10)
        self.axes_upper.grid(True, color='gray')
        self.axes_upper.set_axis_bgcolor('black')
        #self.axes_upper.lines.remove(self.axes_upper_line)
        if self.upperPlot == 'XTraj':
            if not self.get_data(self.sparticle.TRX):
                return
            self.axes_upper.set_xlabel('X', fontsize=9)
            self.axes_upper.set_ylabel('Samples', fontsize=9)
            self.axes_upper.plot(self.ydata, self.tdata)
            #xl0, xl1 = self.sparticle.GetXlim()
            self.axes_upper.set_xlim(self.axes_lower.get_xlim())
        elif self.upperPlot == 'XHist':
            if not self.get_data(self.sparticle.TRX):
                return
            self.axes_upper.set_xlabel('X', fontsize=9)
            self.axes_upper.set_ylabel('X positon Histogram', fontsize=9)
            n, bins, patches = self.axes_upper.hist(self.ydata, 1000, range=self.axes_lower.get_xlim(), histtype='step', normed=False)
            pylab.setp(patches, facecolor='r', alpha=0.75)
        elif self.upperPlot == 'FTraj':
            if not self.get_data(self.sparticle.TR):
                return
            self.axes_upper.set_xlabel('Samples', fontsize=9)
            self.axes_upper.set_ylabel('Potential f(X)', fontsize=9)
            self.axes_upper.plot(self.tdata, self.ydata)
            #yl0, yl1 = self.sparticle.GetYlim()
            self.axes_upper.set_ylim(self.axes_lower.get_ylim())
        elif self.upperPlot == 'FHist':
            if not self.get_data(self.sparticle.TR):
                return
            self.axes_upper.set_xlabel('f(X) Histogram', fontsize=9)
            self.axes_upper.set_ylabel('Potential f(X)', fontsize=9)
            n, bins, patches = self.axes_upper.hist(self.ydata, 1000, range=self.sparticle.GetYlim(), histtype='step', orientation='horizontal' , normed=False)
            pylab.setp(patches, facecolor='r', alpha=0.75)

        pylab.setp(self.axes_upper.get_xticklabels(), fontsize=9)
        pylab.setp(self.axes_upper.get_yticklabels(), fontsize=9)
        #plotting configuration  
        #self.im.set_array(self.sparticle.conf)
        self.axes_lower.lines.remove(self.curXYpoint)
        self.curXYpoint, = self.axes_lower.plot([self.sparticle.dX], [self.sparticle.GetArg(self.sparticle.AI_HAMI)], 'ro')
  
        self.canvas.draw()  
        self.record_figure()
    def on_analysis(self, event):
        if  len(self.ydata) == 0:
            return

        if self.upperPlot == 'XTraj' or self.upperPlot == 'XHist':
            self.draw_ananlysis(r'$\bar X=%g\pm%g$, $\tau_i=%g$', Y=False)
        elif self.upperPlot == 'FTraj' or self.upperPlot == 'FHist':
            self.draw_ananlysis(r'$\bar F=%g\pm%g$, $\tau_i=%g$', Y=True)
        return
    def draw_ananlysis(self, labels, Y):
            m, s, merr , ptau, taui, tauCorr = PyJackknife_MeanStdErrTau_D(self.ydata, 10, 80)
            if taui == -1:
                return 
            print "statistic of m, s, merr , ptau, taui, tauCorr:", m, s, merr , ptau, taui
            print tauCorr
            if Y:
                rect = matplotlib.patches.Rectangle((0, m - merr), width=1, height=2 * merr, transform=self.transY, facecolor='yellow', alpha=0.5)
                l1 = matplotlib.lines.Line2D([0, 1], [m, m], transform=self.transY, linestyle="-.", color="red")
            else:
                rect = matplotlib.patches.Rectangle((m - merr, 0), width=2 * merr, height=1, transform=self.transX, facecolor='yellow', alpha=0.5)
                l1 = matplotlib.lines.Line2D([m, m], [0, 1], transform=self.transX, linestyle="-.", color="red")
  
            self.axes_upper.text(0.05, 0.9, labels % (m, merr, ptau), horizontalalignment='left', verticalalignment='top', color="red", size=10, transform=self.axes_upper.transAxes)
            self.axes_upper.lines.extend([l1])
            self.axes_upper.add_patch(rect)
            self.canvas.draw()   
            self.record_figure()
    def on_potential(self, event):
        self.curchoice = self.radio_box_potential.GetSelection()
        if self.curchoice == 0:
            self.sparticle.Potential = "DB"
        elif self.curchoice == 1:
            self.sparticle.Potential = "Mul"
        else:
            file_choices = "Potential function (*.py)|*.py"
            self.ConfPath = ""
            dlg = wx.FileDialog(
                self,
                message="Open Potential function file...",
                defaultDir=getcwd(),
                defaultFile="",
                wildcard=file_choices,
                style=wx.OPEN)
        
            if dlg.ShowModal() == wx.ID_OK:
                self.ConfPath = dlg.GetPath()
                self.sparticle.Potential = "Custom"
        self.on_reset(event)
        #event.Skip()
  
    def on_upperplot(self, event):
        curchoice = self.radio_box_upperplot.GetSelection()
        if curchoice == 0:
            self.upperPlot = "XTraj"
        elif curchoice == 1:
            self.upperPlot = "XHist"
        elif curchoice == 2:
            self.upperPlot = "FTraj"
        elif curchoice == 3:
            self.upperPlot = "FHist"
        #event.Skip()
    def on_range(self, event):
        self.range = split('[,;: ]?', self.text_ctrl_range.GetValue())
  
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
        self.sparticle.RefreshSamp = f;    
        
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
        self.sparticle.Stride = f;    
    def on_temperature(self, event):
        try:
            f = float(self.text_ctrl_temperature.GetValue())
        except:
            self.mainFrame.refresh_status_message("Unknown temperature:", self.text_ctrl_temperature.GetValue())
            return
        if f < 0:
            self.mainFrame.refresh_status_message("Too small value %d" % f + " and reset to 0.1")
            f = 0.1
            return
        self.sparticle.SetArg(self.sparticle.AI_TEMPER, f)
    def on_wolff_Eavg(self, event):
        try:
            f = int(self.text_ctrl_wolff_Eavg.GetValue())
        except:
            self.mainFrame.refresh_status_message("Unknown wolff_Eavg:", self.text_ctrl_wolff_Eavg.GetValue())
            return
        if f < 2:
            self.mainFrame.refresh_status_message("Too small value %d" % f + " and reset to 2")
            f = 2
            return
        self.sparticle.nGCEWF_AvgC = f;      
        
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
        
        file_choices = "traj (*.txt)|*.txt"
        if self.sparticle.Potential == "DB":
            trajfile = "nvtspak%lg_%lg~%lg_%ld_%dTraj.txt" % (self.sparticle.dX, self.sparticle.dOffset,
                self.sparticle.GetArg(self.sparticle.AI_TEMPER), len(self.sparticle.TR), self.sparticle.Stride);
        elif self.sparticle.Potential == "Mul":
            trajfile = "nvtsmp%lg_%lg~%lg_%ld_%dTraj.txt" % (self.sparticle.dX, self.sparticle.dOffset,
                self.sparticle.GetArg(self.sparticle.AI_TEMPER), len(self.sparticle.TR), self.sparticle.Stride);
                        
        dlg = wx.FileDialog(self, message="Save trajectory file...",
                defaultDir=getcwd(), defaultFile=trajfile, wildcard=file_choices,
                style=wx.SAVE)
                
        if dlg.ShowModal() == wx.ID_OK:
                trajfile = dlg.GetPath()
                f = open(trajfile, 'w')
                if f:
                    for n in range(len(self.sparticle.TR)):
                        f.write(str(self.sparticle.TRX[n]) + "\t" + str(self.sparticle.TR[n]) + "\n")

    def save_conf(self):
        if self.keeping: 
            self.message_dlg("Error", "Stop runing first!")
            return
        if not self.reset: 
            self.message_dlg("Error", "No resetting!")
            return 
        
    def record_figure(self):
        if self.record_figure_checked:
            figfile = "spart%lg_%lg~%lg_%ld_%d_%ldplot%s.png" % (self.sparticle.dX, self.sparticle.dOffset, \
            self.sparticle.GetArg(self.sparticle.AI_TEMPER), len(self.sparticle.TR), self.sparticle.Stride, self.sparticle.uCurMove, strftime("%Y-%m-%d@%H~%M~%S"));
            self.fig.savefig(figfile) 
            
            print "save figure to:", figfile
        else:
            pass

