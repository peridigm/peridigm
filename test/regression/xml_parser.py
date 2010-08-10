#! /usr/bin/env python

import string
import xml.parsers.expat

class XML_Parser:
    def __init__(self, *args):
        self.parser = xml.parsers.expat.ParserCreate()
        self.parser.StartElementHandler = self.start_element
        self.parser.EndElementHandler = self.end_element
        self.parser.CharacterDataHandler = self.char_data
        self.parser.buffer_text = True
        self.parser.returns_unicode = False

    def Parse(self, xml_file_name):
        xml_file = open(xml_file_name)
        self.parser.ParseFile(xml_file)
        xml_file.close()

    def start_element(self, name, attrs):
        return

    def end_element(self, name):
        return

    def char_data(self, data):
        return

class Tolerance_XML_Parser(XML_Parser):

    def __init__(self, *args):
        XML_Parser.__init__(self, *args)
        self.tolerances = {}

    def start_element(self, name, attrs):
        if name == "Tolerance":
            self.tolerances[attrs["variable"]] = float(attrs["value"])

class PVTU_XML_Parser(XML_Parser):

    def __init__(self, *args):
        XML_Parser.__init__(self, *args)
        self.current_label = "None"
        self.pieces = []
        self.point_data_names = []

    def start_element(self, name, attrs):
        # Piece
        if name == "Piece":
            self.pieces.append(str(attrs["Source"]))
        # PPointData
        if name == "PPointData":
            self.current_label = name
        # PDataArray within PPointData
        if self.current_label == "PPointData":
            if name == "PDataArray":
                self.point_data_names.append(str(attrs["Name"]))

    def end_element(self, name):
        if name == "PPointData":
            self.current_label = "None"

class VTU_XML_Info_Parser(XML_Parser):

    def __init__(self, *args):
        XML_Parser.__init__(self, *args)
        self.current_label = "None"
        self.num_points = 0
        self.point_data_names = []

    def start_element(self, name, attrs):
        # Piece
        if name == "Piece":
            num_points = attrs['NumberOfPoints']
        # PointData
        if name == "PointData":
            self.current_label = name
        # DataArray within PointData
        if self.current_label == "PointData":
            if name == "DataArray":
                self.point_data_names.append(str(attrs["Name"]))

    def end_element(self, name):
        if name == "PointData":
            self.current_label = "None"

class VTU_XML_Data_Parser(XML_Parser):

    def __init__(self, *args):
        XML_Parser.__init__(self, *args)
        self.data_arrays_to_record = args[0][:]
        self.data_arrays_to_record.append('ID')
        self.current_label = "None"
        self.current_data = "None"
        self.current_data_type = "None"
        self.num_points = 0
        self.point_data = {}
        self.char_buffer = ""

    def start_element(self, name, attrs):
        # Piece
        if name == "Piece":
            self.num_points = int(attrs['NumberOfPoints'])
        # PointData
        if name == "PointData":
            self.current_label = name
        # DataArray within PointData
        if self.current_label == "PointData":
            if name == "DataArray":
                # is the data name one that we want to record?
                if attrs['Name'] in self.data_arrays_to_record:
                    self.current_data = str(attrs['Name'])
                    self.current_data_type = str(attrs['type'])
                    if self.current_data not in self.point_data.keys():
                        self.point_data[self.current_data] = []

    def end_element(self, name):
        if name == "PointData":
            self.current_label = "None"
        if name == "DataArray":
            # process tail end of last data buffer
            last_piece = string.strip(self.char_buffer)
            if len(last_piece) > 0:
                last_piece = float(self.char_buffer)
                if self.current_data_type == "Int64":
                    last_piece = int(self.char_buffer)
                self.point_data[self.current_data].append(last_piece)
            self.current_data = "None"
            self.current_data_type = "None"
            self.char_buffer = ""

    def char_data(self, data):
        # Data is delivered in chunks that can
        # cut ascii representations of floating
        # point numbers into two meaningless pieces

        # So watch out!

        if len(data) == 0:
            return

        # record the data
        if self.current_label == "PointData" and self.current_data != "None":
            # to avoid chopped off numbers, store tail end of data
            # and tack it on to the next stream
            my_data = self.char_buffer + data
            vals = string.splitfields(my_data)
            if len(vals) == 0:
                # data was not empty but contained only whitespace
                self.char_buffer = ' '
                return
            last_whitespace = string.rfind(my_data, ' ')
            last_newline = string.rfind(my_data, 'n')
            if last_newline > last_whitespace:
                last_whitespace = last_newline
            self.char_buffer = my_data[last_whitespace:]
            if len(string.strip(self.char_buffer)) > 0:
                vals = vals[:-1]

            data_ptr = self.point_data[self.current_data]
            if self.current_data_type == "Int64":
                data_ptr.extend(map(int, vals))
            elif self.current_data_type == "Float64":
                data_ptr.extend(map(float, vals))
            else:
                print "\nxml_parser.py INVALID DATA TYPE", self.current_data_type, "\n"
                sys.exit(1)

        return

class VTU_XML_Points_Parser(XML_Parser):

    def __init__(self, *args):
        XML_Parser.__init__(self, *args)
        self.current_label = "None"
        self.current_data = "None"
        self.current_data_type = "None"
        self.num_points = 0
        self.points = []
        self.char_buffer = ""

    def start_element(self, name, attrs):
        # Piece
        if name == "Piece":
            self.num_points = int(attrs['NumberOfPoints'])
        # Points
        if name == "Points":
            self.current_label = name
        # DataArray within Points
        if self.current_label == "Points":
            if name == "DataArray":
                # check to make sure Name="Points"
                if attrs['Name'] == "Points":
                    self.current_data = str(attrs['Name'])
                    self.current_data_type = str(attrs['type'])

    def end_element(self, name):
        if name == "Points":
            self.current_label = "None"
        if name == "DataArray":
            # process tail end of last data buffer
            last_piece = string.strip(self.char_buffer)
            if len(last_piece) > 0:
                last_piece = float(self.char_buffer)
                self.points.append(last_piece)
            self.current_data = "None"
            self.current_data_type = "None"
            self.char_buffer = ""

    def char_data(self, data):
        # Data is delivered in chunks that can
        # cut ascii representations of floating
        # point numbers into two meaningless pieces

        # So watch out!

        if len(data) == 0:
            return

        # record the data
        if self.current_label == "Points" and self.current_data != "None":
            # to avoid chopped off numbers, store tail end of data
            # and tack it on to the next stream
            my_data = self.char_buffer + data
            vals = string.splitfields(my_data)
            if len(vals) == 0:
                # data was not empty but contained only whitespace
                self.char_buffer = ' '
                return
            last_whitespace = string.rfind(my_data, ' ')
            last_newline = string.rfind(my_data, 'n')
            if last_newline > last_whitespace:
                last_whitespace = last_newline
            self.char_buffer = my_data[last_whitespace:]
            if len(string.strip(self.char_buffer)) > 0:
                vals = vals[:-1]

            # convert the strings to floats and append them to the points array
            self.points.extend(map(float, vals))

        return
