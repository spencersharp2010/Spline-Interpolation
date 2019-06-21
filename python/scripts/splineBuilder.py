import pysplinekernel
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# CLASS DEFINITION
# -----------------------------------------------------------------------------
class SplineBuilder:
    # MEMBER FUNCTIONS --------------------------------------------------------
    def __init__(self):
        #
        self.fig = plt.figure()
        self.ax  = self.fig.add_subplot(111)
        self.title = ' ESC      : reset \n 1..9      : set polynomial degree \n SPACE  : toggle editor mode \n ENTER  : toggle control polygon'
        #
        self.reset()
        
    def reset(self):
        self.ax.cla()
        # Spline setup
        self.p = 3
        self.interpolationPoints = [[],[]]
        # Figure setup
        self.ax.set_title(self.title, loc='left', fontsize=14)
        self.fig.tight_layout()
        # Resize plot
        self.ax.set_xlim( left=0, right=1, auto=False )
        self.ax.set_ylim( bottom=0, top=1, auto=False )
        # Figure objects
        self.spline, = self.ax.plot(0,0,' ')
        self.polygon, = self.ax.plot(0,0,' ')
        self.cursor, = self.ax.plot(0,0,' ')
        # Initialize settings
        self.showControlPolygon = False
        self.inEditor = True
        
    def drawCursor(self, xCursor, yCursor):
        # Draw point at cursor location
        if len(xCursor) is not 0:
            self.cursor.remove()
            self.cursor, = self.ax.plot( xCursor[0], yCursor[0], 'b+' )
        
    def drawSpline(self, xCursor, yCursor):
        # Track cursor location
        self.drawCursor(xCursor,yCursor)
        # Check if number of points is sufficient
        if len(self.interpolationPoints[0]) + len(xCursor)<self.p+1:
            return
        # Get spline
        controlPoints, knotVector = pysplinekernel.interpolateWithBSplineCurve( [self.interpolationPoints[0] + xCursor, self.interpolationPoints[1] + yCursor], self.p )
        # Get parameters (10 per segments)
        t = np.linspace(0.0, 1.0, (len(self.interpolationPoints[0])-1)*10 )
        # Remove previous spline and draw new one
        xc, yc = pysplinekernel.evaluate2DCurve( t, controlPoints[0], controlPoints[1], knotVector )
        self.spline.remove()
        self.spline, = self.ax.plot( xc, yc, 'b' )
        # Draw control polygon
        self.drawControlPolygon(controlPoints)
        #Update figure
        self.fig.canvas.draw()
    
    def drawControlPolygon(self, controlPoints):
        if self.showControlPolygon:
            self.polygon.remove()
            self.polygon, = self.ax.plot(controlPoints[0], controlPoints[1],'r')
        
    # EVENT HANDLERS ----------------------------------------------------------
    def onClick(self,event):
        if not self.inEditor:
            return
        # check if within figure
        if event.inaxes != self.spline.axes: 
            return
        # Check if click location on previous point
        if len(self.interpolationPoints[0]) != 0:
            if event.xdata==self.interpolationPoints[0][-1] and event.ydata==self.interpolationPoints[1][-1]:
                return
        # Append interpolation points
        self.interpolationPoints[0].append(event.xdata)
        self.interpolationPoints[1].append(event.ydata)
        # Draw new point
        self.ax.plot(event.xdata, event.ydata, 'bo')
        self.fig.canvas.draw()
        
    def onMotion(self,event):
        if not self.inEditor:
            return
        # check if within figure
        if event.inaxes != self.spline.axes: 
            return
        # Check if cursor location on previous point
        if len(self.interpolationPoints[0]) != 0:
            if event.xdata==self.interpolationPoints[0][-1] and event.ydata==self.interpolationPoints[1][-1]:
                return
        # Draw spline
        self.drawSpline([event.xdata], [event.ydata])
        self.fig.canvas.draw()
        
    def onKeyPress(self, event):
        # ESC - reset
        if event.key == 'escape':
            self.reset()
        # ENTER - show/hide control polygon
        elif event.key == 'enter':
            if self.showControlPolygon:
                self.polygon.remove()
                self.polygon, = self.ax.plot(0,0,' ')
            self.showControlPolygon = not self.showControlPolygon
            if self.inEditor:
                self.drawSpline([event.xdata],[event.ydata])
            else:
                self.drawSpline([],[])
        # SPACE - stop/start drawing
        elif event.key == ' ':
            self.inEditor = not self.inEditor
            if self.inEditor:
                # Delete temporary interpolation points
                self.interpolationPoints[0] = self.interpolationPoints[0][:-1]
                self.interpolationPoints[1] = self.interpolationPoints[1][:-1]
                # Draw
                self.drawSpline([event.xdata],[event.ydata])
            else:
                # Temporarily append interpolation points
                self.interpolationPoints[0].append(event.xdata)
                self.interpolationPoints[1].append(event.ydata)
                # Draw
                self.drawSpline([],[])
        # NUMBERS - polynomial degree
        elif 49 <= ord(event.key) and ord(event.key) <= 57:
            self.p = int(event.key)
            if self.inEditor:
                self.drawSpline([event.xdata],[event.ydata])
            else:
                self.drawSpline([],[])
        
        
# -----------------------------------------------------------------------------
# END OF CLASS DEFINITION
# -----------------------------------------------------------------------------

# Instantiate SplineBuilder
builder = SplineBuilder()
    
# Define events
motionEvent = builder.ax.figure.canvas.mpl_connect('motion_notify_event', builder.onMotion)
clickEvent = builder.ax.figure.canvas.mpl_connect('button_press_event', builder.onClick)
keyPressEvent = builder.ax.figure.canvas.mpl_connect('key_press_event', builder.onKeyPress)

plt.show()