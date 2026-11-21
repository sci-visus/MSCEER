#!/usr/bin/env python

"""
"""

# tk not integrated with vtk natively
import vtk
import numpy as np
import os
import cv2
from PIL import Image, ImageOps
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import math
import warnings
import tkinter as tk
from tkinter import filedialog
import numpy as np
import random

from tkinter import ttk
from PIL import Image, ImageTk
import time

input_file_name = "C:/Users/jediati/Desktop/JEDIATI/builds/test_GradIntegrator/extract2dridgegraph/Release/optic1.tif_invrt_smoothed_565x584.raw"

#============= GRAPH HERE ==================

class GraphNode:
    def __init__(self, id, pos) :
        self.id = id
        self.pos = pos
        self.arcs = []

    def AddArc(self, arc):
        self.arcs.append(arc)

class GraphArc:

    def __init__(self, node1, node2, polyline):
        self.nodes = [node1, node2]
        self.polyline = polyline
        node1.AddArc(self)
        node2.AddArc(self)


class PolyLineGraph :
    def __init__(self):
        self.nodes = {}
        self.arcs = []

    def AddNode(self, id, pos):
        self.nodes[id] = GraphNode(id, pos)

    def AddArcById(self, node1id, node2id, polyline):
        self.arcs.append(GraphArc(self.nodes[node1id], self.nodes[node2id], polyline))

    def GetNode(self, id):
        return self.nodes[id]

    def GetArcs(self):
        return self.arcs

    def GetNodes(self):
        return self.nodes

    def LoadFromBaseName(self, filenamebase):
        node_name = filenamebase + '.nodes.txt'
        arcs_name = filenamebase + '.arcs.txt'

        file_nodes = open(node_name, 'r')
        # read all the lines
        Lines = file_nodes.readlines()
        file_nodes.close()

        # fill the arrays from the lines read in
        for line in Lines:
            tokens = line.split(',')
            self.AddNode(int(tokens[0]), [float(tokens[1]), float(tokens[2])])

        file_arcs = open(arcs_name, 'r')
        # read all the lines
        Lines = file_arcs.readlines()
        file_arcs.close()

        # fill the arrays from the lines read in
        for line in Lines:
            tokens = line.split(',')
            positions = [float(x) + 0.5 for x in tokens[3:]]
            pos_array = np.array(positions, dtype=np.float).reshape((-1, 2))
            node1id = int(tokens[1])
            node2id = int(tokens[2])
            self.AddArcById(node1id, node2id, pos_array)

#========== END GRAPH HERE ===========

#========== INTERFACE HERE ============
class AutoScrollbar(ttk.Scrollbar):
    """ A scrollbar that hides itself if it's not needed. Works only for grid geometry manager """
    def set(self, lo, hi):
        if float(lo) <= 0.0 and float(hi) >= 1.0:
            self.grid_remove()
        else:
            self.grid()
            ttk.Scrollbar.set(self, lo, hi)

    def pack(self, **kw):
        raise tk.TclError('Cannot use pack with the widget ' + self.__class__.__name__)

    def place(self, **kw):
        raise tk.TclError('Cannot use place with the widget ' + self.__class__.__name__)

class Zoom_Advanced(ttk.Frame):

    def update_image_array(self, numpy_buffer):
        t0 = time.clock()
        self.numpy_buffer = numpy_buffer
        self.show_image()

        print("show(s): ", time.clock() - t0)
    ''' Advanced zoom of the image '''
    def __init__(self, mainframe, numpy_buffer):
        ''' Initialize the main Frame '''
        ttk.Frame.__init__(self, master=mainframe)
        self.master.title('Zoom with mouse wheel')
        # Vertical and horizontal scrollbars for canvas
        vbar = AutoScrollbar(self.master, orient='vertical')
        hbar = AutoScrollbar(self.master, orient='horizontal')
        vbar.grid(row=0, column=1, sticky='ns')
        hbar.grid(row=1, column=0, sticky='we')
        # Create canvas and put image on it
        self.canvas = tk.Canvas(self.master, highlightthickness=0,
                                xscrollcommand=hbar.set, yscrollcommand=vbar.set)
        self.canvas.grid(row=0, column=0, sticky='nswe')
        self.canvas.update()  # wait till canvas is created
        vbar.configure(command=self.scroll_y)  # bind scrollbars to the canvas
        hbar.configure(command=self.scroll_x)
        # Make the canvas expandable
        self.master.rowconfigure(0, weight=1)
        self.master.columnconfigure(0, weight=1)
        # Bind events to the Canvas
        self.canvas.bind('<Configure>', self.show_image)  # canvas is resized
        self.canvas.bind('<ButtonPress-1>', self.move_from)
        self.canvas.bind('<B1-Motion>',     self.move_to)
        self.canvas.bind('<MouseWheel>', self.wheel)  # with Windows and MacOS, but not Linux
        self.canvas.bind('<Button-5>',   self.wheel)  # only with Linux, wheel scroll down
        self.canvas.bind('<Button-4>',   self.wheel)  # only with Linux, wheel scroll up
        self.numpy_buffer = numpy_buffer
        self.image = ImageOps.flip(Image.fromarray(self.numpy_buffer))  # open image

        self.width, self.height = self.image.size
        self.imwidth, self.imheight = self.image.size
        self.imscale = 1.0  # scale for the canvaas image
        self.delta = 1.3  # zoom magnitude
        # Put image into container rectangle and use it to set proper coordinates to the image
        self.container = self.canvas.create_rectangle(0, 0, self.width, self.height, width=0)
        # Plot some optional random rectangles for the test purposes
        # minsize, maxsize, number = 5, 20, 10
        # for n in range(number):
        #     x0 = random.randint(0, self.width - maxsize)
        #     y0 = random.randint(0, self.height - maxsize)
        #     x1 = x0 + random.randint(minsize, maxsize)
        #     y1 = y0 + random.randint(minsize, maxsize)
        #     color = ('red', 'orange', 'yellow', 'green', 'blue')[random.randint(0, 4)]
        #     self.canvas.create_rectangle(x0, y0, x1, y1, fill=color, activefill='black')
        self.show_image()

    def scroll_y(self, *args, **kwargs):
        ''' Scroll canvas vertically and redraw the image '''
        self.canvas.yview(*args, **kwargs)  # scroll vertically
        self.show_image()  # redraw the image

    def scroll_x(self, *args, **kwargs):
        ''' Scroll canvas horizontally and redraw the image '''
        self.canvas.xview(*args, **kwargs)  # scroll horizontally
        self.show_image()  # redraw the image

    def move_from(self, event):
        ''' Remember previous coordinates for scrolling with the mouse '''
        self.canvas.scan_mark(event.x, event.y)

    def move_to(self, event):
        ''' Drag (move) canvas to the new position '''
        self.canvas.scan_dragto(event.x, event.y, gain=1)
        self.show_image()  # redraw the image

    def wheel(self, event):
        ''' Zoom with mouse wheel '''
        x = self.canvas.canvasx(event.x)
        y = self.canvas.canvasy(event.y)
        bbox = self.canvas.bbox(self.container)  # get image area
        if bbox[0] < x < bbox[2] and bbox[1] < y < bbox[3]: pass  # Ok! Inside the image
        else: return  # zoom only inside image area
        scale = 1.0
        # Respond to Linux (event.num) or Windows (event.delta) wheel event
        if event.num == 5 or event.delta == -120:  # scroll down
            i = min(self.width, self.height)
            if int(i * self.imscale) < 30: return  # image is less than 30 pixels
            self.imscale /= self.delta
            scale        /= self.delta
        if event.num == 4 or event.delta == 120:  # scroll up
            i = min(self.canvas.winfo_width(), self.canvas.winfo_height())
            if i < self.imscale: return  # 1 pixel is bigger than the visible area
            self.imscale *= self.delta
            scale        *= self.delta
        self.canvas.scale('all', x, y, scale, scale)  # rescale all canvas objects
        self.show_image()

    def show_image(self, event=None):
        ''' Show image on the Canvas '''
        bbox1 = self.canvas.bbox(self.container)  # get image area
        # Remove 1 pixel shift at the sides of the bbox1
        bbox1 = (bbox1[0] + 1, bbox1[1] + 1, bbox1[2] - 1, bbox1[3] - 1)
        bbox2 = (self.canvas.canvasx(0),  # get visible area of the canvas
                 self.canvas.canvasy(0),
                 self.canvas.canvasx(self.canvas.winfo_width()),
                 self.canvas.canvasy(self.canvas.winfo_height()))
        bbox = [min(bbox1[0], bbox2[0]), min(bbox1[1], bbox2[1]),  # get scroll region box
                max(bbox1[2], bbox2[2]), max(bbox1[3], bbox2[3])]
        if bbox[0] == bbox2[0] and bbox[2] == bbox2[2]:  # whole image in the visible area
            bbox[0] = bbox1[0]
            bbox[2] = bbox1[2]
        if bbox[1] == bbox2[1] and bbox[3] == bbox2[3]:  # whole image in the visible area
            bbox[1] = bbox1[1]
            bbox[3] = bbox1[3]
        self.canvas.configure(scrollregion=bbox)  # set scroll region
        x1 = max(bbox2[0] - bbox1[0], 0)  # get coordinates (x1,y1,x2,y2) of the image tile
        y1 = max(bbox2[1] - bbox1[1], 0)
        x2 = min(bbox2[2], bbox1[2]) - bbox1[0]
        y2 = min(bbox2[3], bbox1[3]) - bbox1[1]
        if int(x2 - x1) > 0 and int(y2 - y1) > 0:  # show image if it in the visible area
            x = min(int(x2 / self.imscale), self.width)   # sometimes it is larger on 1 pixel...
            y = min(int(y2 / self.imscale), self.height)  # ...and sometimes not
            image = self.image.crop((int(x1 / self.imscale), int(y1 / self.imscale), x, y))
            imagetk = ImageTk.PhotoImage(image.resize((int(x2 - x1), int(y2 - y1))))
            imageid = self.canvas.create_image(max(bbox2[0], bbox1[0]), max(bbox2[1], bbox1[1]),
                                               anchor='nw', image=imagetk)
            self.canvas.lower(imageid)  # set image into background
            self.canvas.imagetk = imagetk  # keep an extra reference to prevent garbage-collection
    """ Display and zoom image """
    def scale_to_box(self, x, y):
        """ Checks if the point (x,y) is outside the image area """
        bbox = self.canvas.coords(self.container)  # get image area
        # parametrize x
        rx = (x - bbox[0]) / float(bbox[2] - bbox[0])
        # parametrize y
        ry = (y - bbox[1]) / float(bbox[3] - bbox[1])
        return rx, ry

    # get the bounding box (x0, y0, x1, y1)
    def get_bbox(self):
        return self.canvas.coords(self.container)

    #get the pixel coordinates of normalized coordinates
    def screen_to_pixel_coordinates(self, x, y):
        nx, ny = self.scale_to_box(x,y)
        return nx * self.imwidth, ny * self.imheight



class ImageOverlayScene(Zoom_Advanced):
    def make_oval(self, nx, ny, r):
        bbox = self.canvas.coords(self.container)
        dx = bbox[2] - bbox[0]
        dy = bbox[3] - bbox[1]
        x0 = (float(nx) / self.imwidth) * dx + bbox[0] - r
        y0 = (float(ny) / self.imheight) * dy + bbox[1] - r
        x1 = (float(nx) / self.imwidth) * dx + bbox[0] + r
        y1 = (float(ny) / self.imheight) * dy + bbox[1] + r

        id = self.canvas.create_oval(x0,y0,x1, y1, tag="Node", fill="White", activefill="Red", activeoutline="Yellow")
        return id


    def __motion(self, event):
        #print("moved", event.x, event.y)
        #closest = self.canvas.find_closest(event.x, event.y)
        #print("closest=", closest)
        x = self.canvas.canvasx(event.x)  # get coordinates of the event on the canvas
        y = self.canvas.canvasy(event.y)
        #if self.outside(x, y): return  # zoom only inside image area
        rx, ry = self.scale_to_box(x,y)
        print(rx * self.imwidth, ry * self.imheight)
        self.mouse_over_x, self.mouse_over_y = rx * self.imwidth, ry * self.imheight

        t0  = time.clock()
        buffer = self.buffer_renderer.Render()
        self.update_image_array(buffer)
        print("motion (s):", time.clock() - t0)

    def __init__(self, placeholder, buff_renderer):
        self.buffer_renderer = buff_renderer
        Zoom_Advanced.__init__(self, placeholder, buff_renderer.Render())
        self.canvas.bind('<Double-Button-1>', self.__shift_click_draw) # to queries
        self.canvas.bind('<Motion>', self.__motion)
        self.__list_of_points = []
        self.__list_of_ovals = []
        self.__radius = 5
        self.mouse_over_x, self.mouse_over_y = 0, 0


    def __shift_click_draw(self, event):
        print("shift click")
        x = self.canvas.canvasx(event.x)  # get coordinates of the event on the canvas
        y = self.canvas.canvasy(event.y)
        #if self.outside(x, y): return  # zoom only inside image area
        rx, ry = self.screen_to_pixel_coordinates(x,y)
        print(" -- in-image:", x, y, "==>", rx, ry)
        self.__list_of_points.append( (rx , ry))
        self.__list_of_ovals.append(self.make_oval(rx, ry, self.__radius))

    def redraw_figures(self):
        bbox = self.canvas.coords(self.container)
        dx = bbox[2] - bbox[0]
        dy = bbox[3] - bbox[1]
        #print(self.canvas.gettags("Node"))
        for oval in self.canvas.find_withtag("Node"):
            x0, y0, x1, y1 = self.canvas.coords(oval)
            print("old=", x0, y0, x1, y1)
            mx = (x0 + x1) / 2
            my = (y0 + y1) / 2

            x0 = mx - self.__radius
            y0 = my - self.__radius
            x1 = mx + self.__radius
            y1 = my + self.__radius

            print("new=", x0, y0, x1, y1)
            self.canvas.coords(oval, x0, y0, x1, y1)

#============== END INTERFACE HERE ============

#============== VTK RENDER TO ARRAY HERE ========

class VTKBufferRenderer:
    def __init__(self):
        self.renderWindow = vtk.vtkRenderWindow()
        self.renderWindow.SetOffScreenRendering(1)

        self.windowToImageFilter = vtk.vtkWindowToImageFilter()
        self.windowToImageFilter.SetInput(self.renderWindow)

        self.renderers = []
        self.writer = vtk.vtkPNGWriter()
        self.writer.SetWriteToMemory(1)
        self.writer.SetInputConnection(self.windowToImageFilter.GetOutputPort())

    def SetSize(self, width, height):
        self.width = width
        self.height = height
        print("setting width, height:", width, height)
        self.renderWindow.SetSize(width , height )

        self.camera = vtk.vtkCamera()
        self.xFocal = (width - 1) / 2.0
        self.yFocal = (height - 1) / 2.0
        self.camera.SetFocalPoint(self.xFocal, self.yFocal, 0)
        self.camera.SetPosition(self.xFocal, self.yFocal, 1)
        self.camera.SetViewUp(0,1,0)
        self.camera.ParallelProjectionOn()
        self.camera.SetParallelScale(max(width, height)/2)
        #self.windowToImageFilter.SetScale(2,2)


    def AddRenderer(self, renderer):
        self.renderWindow.AddRenderer(renderer)
        self.renderers.append(renderer)


    def Render(self):
        t0 = time.clock()

        #camera = self.renderers[0].GetActiveCamera()
        self.renderers[0].SetActiveCamera(self.camera)
        self.renderWindow.Render()
        self.windowToImageFilter.Update()
        test = self.windowToImageFilter.GetOutput()
        #print(test)
        width, height, depth = test.GetDimensions()
        #print("whd: ", width, height, depth)
        #print(self.writer)
        #self.writer.Write()
        #data = bytes(memoryview(self.writer.GetResult()))
        vtk_array = test.GetPointData().GetScalars()
        components = vtk_array.GetNumberOfComponents()
        arr = vtk_to_numpy(vtk_array).reshape(height, width, components)
        print("array shape:", arr.shape)
        print("render(s):", time.clock() - t0)
        return arr

def getDimsFromFilename(filename_raw) :
    dirname, filename = os.path.split(filename_raw)
    print(dirname, filename)
    filenamebase = os.path.splitext(filename)[0]
    dim_string_list = filenamebase.split("_")[-1].split("x")
    return [int(x) for x in dim_string_list]

class VTKImageActor:
    def LoadRawImage(self, raw_filename ):
        self.raw_filename = raw_filename
        self.dims = getDimsFromFilename(raw_filename)
        if len(self.dims) == 2:
            self.dims.append(1)
        self.reader = vtk.vtkImageReader()
        self.reader.SetDataScalarType(vtk.VTK_FLOAT)  # unsigned int8
        self.reader.SetFileName(raw_filename)
        self.reader.SetNumberOfScalarComponents(1)
        self.reader.SetFileDimensionality(3)
        self.reader.SetDataByteOrderToLittleEndian()
        extent_list = []
        for x in self.dims:
            extent_list.append(0)
            extent_list.append(x-1)
        print("data extent: ", extent_list)
        self.reader.SetDataExtent(*extent_list)  #  image size 488*488*332
        self.reader.SetDataSpacing(1,1,1)  # Volume Pixel
        self.reader.Update()
        self.scalarRange =  self.reader.GetOutput().GetPointData().GetScalars().GetRange()
        print("Range:", self.scalarRange)

        # Work with triple images.
        self.cast = vtk.vtkImageCast()
        self.cast.SetInputConnection(self.reader.GetOutputPort())
        self.cast.SetOutputScalarTypeToDouble()
        self.cast.Update()


        self.colorWindow = (self.scalarRange[1] - self.scalarRange[0])
        self.colorLevel = self.colorWindow

        # Map the image through the lookup table.
        self.originalColor = vtk.vtkImageMapToWindowLevelColors()
        self.originalColor.SetWindow(self.colorWindow)
        self.originalColor.SetLevel(self.colorLevel)
        self.originalColor.SetInputConnection(self.reader.GetOutputPort())

        # make an actor for it
        self.originalActor = vtk.vtkImageActor()
        self.originalActor.GetMapper().SetInputConnection(self.originalColor.GetOutputPort())
        self.originalActor.GetProperty().SetInterpolationTypeToNearest()
        self.originalActor.SetDisplayExtent(*self.reader.GetDataExtent())

    def GetActor(self):
        return self.originalActor

    def GetDims(self):
        return self.dims


class VTKGraphActor:
    def LoadFromGraph(self, raw_filename ):
        self.raw_filename = raw_filename
        self.dims = getDimsFromFilename(raw_filename)
        if len(self.dims) == 2:
            self.dims.append(1)
        self.plg = PolyLineGraph()
        self.plg.LoadFromBaseName(raw_filename)

        self.raw_filename = raw_filename

        self.all_points = []
        for arc in self.plg.GetArcs():
            self.all_points.extend(arc.polyline.tolist())
        points = vtk.vtkPoints()
        self.all_points = [[x[0], self.dims[1] - x[1],0.1] for x in self.all_points]
        #print(self.all_points)
        points.SetData(numpy_to_vtk(np.asarray(self.all_points)))
        point_offset = 0
        cells = vtk.vtkCellArray()
        self.all_vtk_polylines = []
        for arc in self.plg.GetArcs():
            polyline = vtk.vtkPolyLine()
            numpoints = arc.polyline.shape[0]
            polyline.GetPointIds().SetNumberOfIds(numpoints)
            print("arc has points: ", numpoints)
            for i in range(numpoints):
                polyline.GetPointIds().SetId(i,i + point_offset)
            cells.InsertNextCell(polyline)
            point_offset = point_offset + numpoints
        self.polyData = vtk.vtkPolyData()
        self.polyData.SetPoints(points)
        self.polyData.SetLines(cells)
        # Setup actor and mapper
        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInputData(self.polyData)

        self.actor = vtk.vtkActor()
        self.actor.SetMapper(self.mapper)
        self.actor.GetProperty().SetColor(1.0,0.5,0.1)
        self.actor.GetProperty().SetLineWidth(1)
        self.actor.ForceTranslucentOn()

    def GetActor(self):
        return self.actor

    def GetDims(self):
        return self.dims

#=========== END VTK RENDER TO IMAGE HERE =============




imren = VTKImageActor()
imren.LoadRawImage(input_file_name)

graph = VTKGraphActor()
graph.LoadFromGraph(input_file_name)

# originalRenderer.SetViewport([0,0,1,1])
originalRenderer = vtk.vtkRenderer()
originalRenderer.AddActor(imren.GetActor())
originalRenderer.AddActor(graph.GetActor())

# originalRenderer.GetActiveCamera().Dolly(-5.5)
#
#originalRenderer.GetActiveCamera().ParallelProjectionOn()


buff_renderer = VTKBufferRenderer()
buff_renderer.SetSize(*(imren.GetDims())[:2])
buff_renderer.AddRenderer(originalRenderer)
buff_renderer.Render()
# originalRenderer.ResetCameraClippingRange()
# camera = originalRenderer.GetActiveCamera()
# camera.ParallelProjectionOn()
# camera.SetParallelScale(2)


#img = Image.fromarray(buffer)
#img.show()
root = tk.Tk()
root.geometry(str(imren.GetDims()[0] + 20) + "x" + str(imren.GetDims()[1]+ 20))
print("making scene")
app = ImageOverlayScene(root, buff_renderer)
print("done making scene")
t0 = time.clock()

root.mainloop()

print("got here")
# print("Got Here")
# gather_arcs = [arc.polyline for arc in plg.GetArcs()]
#
# image = np.zeros((1028,1024,3), np.uint8)
# image = cv2.polylines(image, gather_arcs, False, (255,0,0), 2)
# # Displaying the image
# while (1):
#
#     cv2.imshow('image', image)
#     if cv2.waitKey(20) & 0xFF == 27:
#         break
#
# cv2.destroyAllWindows()


# def main():
#     # colors = vtk.vtkNamedColors()
#
#     fileName = "C:/Users/jediati/Desktop/JEDIATI/docs/figures/2_pari_full_gbtest_small_n_80x187x357-1.png"
#
#     # Read the image.
#     readerFactory = vtk.vtkImageReader2Factory()
#     reader = readerFactory.CreateImageReader2(fileName)
#     reader.SetFileName(fileName)
#     reader.Update()
#
#     scalarRange = [0] * 2
#     scalarRange[0] = reader.GetOutput().GetPointData().GetScalars().GetRange()[0]
#     scalarRange[1] = reader.GetOutput().GetPointData().GetScalars().GetRange()[1]
#     print("Range:", scalarRange)
#     middleSlice = 0
#
#     # Work with triple images.
#     cast = vtk.vtkImageCast()
#     cast.SetInputConnection(reader.GetOutputPort())
#     cast.SetOutputScalarTypeToDouble()
#     cast.Update()
#
#     laplacian = vtk.vtkImageLaplacian()
#     laplacian.SetInputConnection(cast.GetOutputPort())
#     laplacian.SetDimensionality(2)
#
#     enhance = vtk.vtkImageMathematics()
#     enhance.SetInputConnection(0, cast.GetOutputPort())
#     enhance.SetInputConnection(1, laplacian.GetOutputPort())
#     enhance.SetOperationToSubtract()
#
#     colorWindow = (scalarRange[1] - scalarRange[0])
#     colorLevel = colorWindow / 2
#
#     # Map the image through the lookup table.
#     originalColor = vtk.vtkImageMapToWindowLevelColors()
#     originalColor.SetWindow(colorWindow)
#     originalColor.SetLevel(colorLevel)
#     originalColor.SetInputConnection(reader.GetOutputPort())
#
#     originalActor = vtk.vtkImageActor()
#     originalActor.GetMapper().SetInputConnection(originalColor.GetOutputPort())
#     originalActor.GetProperty().SetInterpolationTypeToNearest()
#     originalActor.SetDisplayExtent(
#         reader.GetDataExtent()[0], reader.GetDataExtent()[1],
#         reader.GetDataExtent()[2], reader.GetDataExtent()[3],
#         middleSlice, middleSlice)
#
#     laplacianColor = vtk.vtkImageMapToWindowLevelColors()
#     laplacianColor.SetWindow(1000)
#     laplacianColor.SetLevel(0)
#     laplacianColor.SetInputConnection(laplacian.GetOutputPort())
#
#     laplacianActor = vtk.vtkImageActor()
#     laplacianActor.GetMapper().SetInputConnection(laplacianColor.GetOutputPort())
#     laplacianActor.GetProperty().SetInterpolationTypeToNearest()
#     laplacianActor.SetDisplayExtent(originalActor.GetDisplayExtent())
#
#     enhancedColor = vtk.vtkImageMapToWindowLevelColors()
#     enhancedColor.SetWindow(colorWindow)
#     enhancedColor.SetLevel(colorLevel)
#     enhancedColor.SetInputConnection(enhance.GetOutputPort())
#
#     enhancedActor = vtk.vtkImageActor()
#     enhancedActor.GetMapper().SetInputConnection(enhancedColor.GetOutputPort())
#     enhancedActor.GetProperty().SetInterpolationTypeToNearest()
#     enhancedActor.SetDisplayExtent(originalActor.GetDisplayExtent())
#
#     # Setup the renderers.
#     originalRenderer = vtk.vtkRenderer()
#     originalRenderer.AddActor(originalActor)
#     laplacianRenderer = vtk.vtkRenderer()
#     laplacianRenderer.AddActor(laplacianActor)
#     enhancedRenderer = vtk.vtkRenderer()
#     enhancedRenderer.AddActor(enhancedActor)
#
#     renderers = list()
#     renderers.append(originalRenderer)
#     renderers.append(laplacianRenderer)
#     renderers.append(enhancedRenderer)
#
#     # Setup viewports for the renderers.
#     rendererSize = 400
#     xGridDimensions = 3
#     yGridDimensions = 1
#
#     renderWindow = vtk.vtkRenderWindow()
#     renderWindow.SetSize(rendererSize * xGridDimensions, rendererSize * yGridDimensions)
#     for row in range(0, yGridDimensions):
#         for col in range(xGridDimensions):
#             index = row * xGridDimensions + col
#             # (xmin, ymin, xmax, ymax)
#             viewport = [float(col) / xGridDimensions, float(yGridDimensions - (row + 1)) / yGridDimensions,
#                         float(col + 1) / xGridDimensions, float(yGridDimensions - row) / yGridDimensions]
#             renderers[index].SetViewport(viewport)
#             renderWindow.AddRenderer(renderers[index])
#     renderWindow.SetWindowName('EnhanceEdges')
#
#     renderWindowInteractor = vtk.vtkRenderWindowInteractor()
#     style = vtk.vtkInteractorStyleImage()
#
#     renderWindowInteractor.SetInteractorStyle(style)
#     renderWindowInteractor.SetRenderWindow(renderWindow)
#
#     # Renderers share one camera.
#     renderWindow.Render()
#     renderers[0].GetActiveCamera().Dolly(1.5)
#     renderers[0].ResetCameraClippingRange()
#
#     for r in range(1, len(renderers)):
#         renderers[r].SetActiveCamera(renderers[0].GetActiveCamera())
#     renderWindowInteractor.Initialize()
#     renderWindowInteractor.Start()
#
#
# def get_program_parameters():
#     import argparse
#     description = 'High-pass filters can extract and enhance edges in an image.'
#     epilogue = '''
#     Subtraction of the Laplacian (middle) from the original image (left) results
#      in edge enhancement or a sharpening operation (right).
#     '''
#     parser = argparse.ArgumentParser(description=description, epilog=epilogue,
#                                      formatter_class=argparse.RawDescriptionHelpFormatter)
#     parser.add_argument('filename', help='FullHead.mhd.')
#     args = parser.parse_args()
#     return args.filename
#
#
# if __name__ == '__main__':
#     main()