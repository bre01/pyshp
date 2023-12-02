
from struct import *
import os
NODATA = -10e38
NULL = 0
POINT = 1
POLYLINE = 3
POLYGON = 5
MULTIPOINT = 8
POINTZ = 11
POLYLINEZ = 13
POLYGONZ = 15
MULTIPOINTZ = 18
POINTM = 21
POLYLINEM = 23
POLYGONM = 25
MULTIPOINTM = 28
MULTIPATCH = 31
izip=zip
xrange=range
SHAPETYPE_LOOKUP = {
    0: 'NULL',
    1: 'POINT',
    3: 'POLYLINE',
    5: 'POLYGON',
    8: 'MULTIPOINT',
    11: 'POINTZ',
    13: 'POLYLINEZ',
    15: 'POLYGONZ',
    18: 'MULTIPOINTZ',
    21: 'POINTM',
    23: 'POLYLINEM',
    25: 'POLYGONM',
    28: 'MULTIPOINTM',
    31: 'MULTIPATCH'}

TRIANGLE_STRIP = 0
TRIANGLE_FAN = 1
OUTER_RING = 2
INNER_RING = 3
FIRST_RING = 4
RING = 5

PARTTYPE_LOOKUP = {
    0: 'TRIANGLE_STRIP',
    1: 'TRIANGLE_FAN',
    2: 'OUTER_RING',
    3: 'INNER_RING',
    4: 'FIRST_RING',
    5: 'RING'}



def u(v, encoding='utf-8', encodingErrors='strict'):
        if isinstance(v, bytes):
            # For python 3 decode bytes to str.
            return v.decode(encoding, encodingErrors)
        elif isinstance(v, str):
            # Already str.
            return v
        elif v is None:
            # Since we're dealing with text, interpret None as ""
            return ""
        else:
            # Force string representation.
            return bytes(v).decode(encoding, encodingErrors)

class Reader(object):
    def __init__(self,*args,**kwargs):
       self.shp=None
       self.shx=None
       self.dbf=None
       self.path=args[0]
       self._file_to_close=[];
       self.encoding = kwargs.pop('encoding', 'utf-8')
       self.encodingErrors = kwargs.pop('encodingErrors', 'strict')
       self.shapeName=None
       self.shp = None
       self.shx = None
       self.dbf = None
       self._files_to_close = []
       self.shapeName = "Not specified"
       self._offsets = []
       self.shpLength = None
       self.numRecords = None
       self.numShapes = None
       self.fields = []
       self.__dbfHdrLength = 0
       self.__fieldLookup = {}
       self.encoding = kwargs.pop('encoding', 'utf-8')
       self.encodingErrors = kwargs.pop('encodingErrors', 'strict')
       self.load(self.path);
    def load(self,shapefile=None):
        (shapeName,ext)=os.path.splitext(shapefile);
        self.shapeName=shapeName
        self.load_shp(shapeName)
        self.load_shx(shapeName)
        self.load_dbf(shapeName)
        if not (self.shp or self.dbf):
            raise Exception("Unable to open dbf or shp");
        if self.shp: 
            self.__shpHeader()
        if self.dbf:
            self.__dbfHeader()
        if self.shx:
            self.__shxHeader()
    def __shxHeader(self):
        """Reads the header information from a .shx file."""
        shx = self.shx
        if not shx:
            raise Exception("Shapefile Reader requires a shapefile or file-like object. (no shx file found")
        # File length (16-bit word * 2 = bytes) - header length
        shx.seek(24)
        shxRecordLength = (unpack(">i", shx.read(4))[0] * 2) - 100
        self.numShapes = shxRecordLength // 8

    def load_shp(self,shapefile_name):
        shp_ext="shp"
        self.shp=open("%s.%s" % (shapefile_name,shp_ext),"rb")
        self._file_to_close.append(self.shp);




    def load_shx(self,file_name):
        self.shx=open("%s.%s"%(file_name,"shx"),"rb")
        self._file_to_close.append(self.shx)
    def load_dbf(self,shapefile_name):
        self.dbf=open("%s.%s" % (shapefile_name,"dbf"),"rb")
        self._file_to_close.append(self.dbf)
        
    def __shpHeader(self):
        shp=self.shp
        shp.seek(24)
        #file length big endian
        self.shpLength= unpack(">i",shp.read(4))[0]*2
        #from word to byte  
        
        #shape type little endian
        shp.seek(32)
        self.shapeType=unpack("<i",shp.read(4))[0]
        self.bbox=_Array('d',unpack("<4d",shp.read(32)))
        self.zbox=_Array('d',unpack("<2d",shp.read(16)))
        self.mbox=[]
        for m in _Array('d',unpack("<2d",shp.read(16))):
            if m > NODATA:
                self.mbox.append(m)
            else:
                self.mbox.append(None)
    @property
    def shapeRecords(self):


        return  ShapeRecords(self.iterShapeRecords())
    
    def iterShapeRecords(self):
        for shape,record in zip(self.iterShapes(),self.iterRecords()):
            yield ShapeRecord(shape=shape,record=record)
    def __dbfHeader(self):
        """Reads a dbf header. Xbase-related code borrows heavily from ActiveState Python Cookbook Recipe 362715 by Raymond Hettinger"""
        if not self.dbf:
            raise Exception("Shapefile Reader requires a shapefile or file-like object. (no dbf file found)")
        dbf = self.dbf
        # read relevant header parts
        dbf.seek(0)
        self.numRecords, self.__dbfHdrLength, self.__recordLength = \
                unpack("<xxxxLHH20x", dbf.read(32))
        # read fields
        numFields = (self.__dbfHdrLength - 33) // 32
        for field in range(numFields):
            fieldDesc = list(unpack("<11sc4xBB14x", dbf.read(32)))
            name = 0
            idx = 0
            if b"\x00" in fieldDesc[name]:
                idx = fieldDesc[name].index(b"\x00")
            else:
                idx = len(fieldDesc[name]) - 1
            fieldDesc[name] = fieldDesc[name][:idx]
            fieldDesc[name] = u(fieldDesc[name], self.encoding, self.encodingErrors)
            fieldDesc[name] = fieldDesc[name].lstrip()
            fieldDesc[1] = u(fieldDesc[1], 'ascii')
            self.fields.append(fieldDesc)
        terminator = dbf.read(1)
        if terminator != b"\r":
            raise Exception("Shapefile dbf header lacks expected terminator. (likely corrupt?)")
        
        # insert deletion field at start
        self.fields.insert(0, ('DeletionFlag', 'C', 1, 0))

        # store all field positions for easy lookups
        # note: fieldLookup gives the index position of a field inside Reader.fields
        self.__fieldLookup = dict((f[0], i) for i, f in enumerate(self.fields))

        # by default, read all fields except the deletion flag, hence "[1:]"
        # note: recLookup gives the index position of a field inside a _Record list
        fieldnames = [f[0] for f in self.fields[1:]]
        fieldTuples,recLookup,recStruct = self.__recordFields(fieldnames)
        self.__fullRecStruct = recStruct
        self.__fullRecLookup = recLookup
    
    def __getFileObj(self,f):
        if self.shp and self.shpLength is None:
            self.load()
        if self.dbf and len(self.fields)  ==0:
            self.load()
        return f
    def __recordFields(self, fields=None):
        """Returns the necessary info required to unpack a record's fields,
        restricted to a subset of fieldnames 'fields' if specified. 
        Returns a list of field info tuples, a name-index lookup dict, 
        and a Struct instance for unpacking these fields. Note that DeletionFlag
        is not a valid field. 
        """
        if fields is not None:
            # restrict info to the specified fields
            # first ignore repeated field names (order doesn't matter)
            fields = list(set(fields))
            # get the struct
            fmt, fmtSize = self.__recordFmt(fields=fields)
            recStruct = Struct(fmt)
            # make sure the given fieldnames exist
            for name in fields:
                if name not in self.__fieldLookup or name == 'DeletionFlag':
                    raise ValueError('"{}" is not a valid field name'.format(name))
            # fetch relevant field info tuples
            fieldTuples = []
            for fieldinfo in self.fields[1:]:
                name = fieldinfo[0]
                if name in fields:
                    fieldTuples.append(fieldinfo)
            # store the field positions
            recLookup = dict((f[0], i) for i,f in enumerate(fieldTuples))
        else:
            # use all the dbf fields
            fieldTuples = self.fields[1:] # sans deletion flag
            recStruct = self.__fullRecStruct
            recLookup = self.__fullRecLookup
        return fieldTuples, recLookup, recStruct
    def __recordFmt(self, fields=None):
        """Calculates the format and size of a .dbf record. Optional 'fields' arg 
        specifies which fieldnames to unpack and which to ignore. Note that this
        always includes the DeletionFlag at index 0, regardless of the 'fields' arg. 
        """
        if self.numRecords is None:
            self.__dbfHeader()
        structcodes = ['%ds' % fieldinfo[2]
                        for fieldinfo in self.fields]
        if fields is not None:
            # only unpack specified fields, ignore others using padbytes (x)
            structcodes = [code if fieldinfo[0] in fields 
                            or fieldinfo[0] == 'DeletionFlag' # always unpack delflag
                            else '%dx' % fieldinfo[2]
                            for fieldinfo,code in zip(self.fields, structcodes)]
        fmt = ''.join(structcodes)
        fmtSize = calcsize(fmt)
        # total size of fields should add up to recordlength from the header
        while fmtSize < self.__recordLength:
            # if not, pad byte until reaches recordlength
            fmt += "x"
            fmtSize += 1
        return (fmt, fmtSize)    


    def iterShapes(self, bbox=None):
        """Returns a generator of shapes in a shapefile. Useful
        for handling large shapefiles.
        To only read shapes within a given spatial region, specify the 'bbox'
        arg as a list or tuple of xmin,ymin,xmax,ymax. 
        """
        shp = self.__getFileObj(self.shp)
        # Found shapefiles which report incorrect
        # shp file length in the header. Can't trust
        # that so we seek to the end of the file
        # and figure it out.
        shp.seek(0,2)
        shpLength = shp.tell()
        shp.seek(100)

        if self.numShapes:
            # Iterate exactly the number of shapes from shx header
            for i in range(self.numShapes):
                # MAYBE: check if more left of file or exit early? 
                shape = self.__shape(oid=i, bbox=bbox)
                if shape:
                    yield shape
        else:
            # No shx file, unknown nr of shapes
            # Instead iterate until reach end of file
            # Collect the offset indices during iteration
            i = 0
            offsets = []
            pos = shp.tell()
            while pos < shpLength:
                offsets.append(pos)
                shape = self.__shape(oid=i, bbox=bbox)
                pos = shp.tell()
                if shape:
                    yield shape
                i += 1
            # Entire shp file consumed
            # Update the number of shapes and list of offsets
            assert i == len(offsets)
            self.numShapes = i 
            self._offsets = offsets
    def __shape(self, oid=None, bbox=None):
        """Returns the header info and geometry for a single shape."""
        f = self.__getFileObj(self.shp)
        record = Shape(oid=oid)
        nParts = nPoints = zmin = zmax = mmin = mmax = None
        (recNum, recLength) = unpack(">2i", f.read(8))
        # Determine the start of the next record
        next = f.tell() + (2 * recLength)
        shapeType = unpack("<i", f.read(4))[0]
        record.shapeType = shapeType
        # For Null shapes create an empty points list for consistency
        if shapeType == 0:
            record.points = []
        # All shape types capable of having a bounding box
        elif shapeType in (3,5,8,13,15,18,23,25,28,31):
            record.bbox = _Array('d', unpack("<4d", f.read(32)))
            # if bbox specified and no overlap, skip this shape
            if bbox is not None and not bbox_overlap(bbox, record.bbox):
                # because we stop parsing this shape, skip to beginning of
                # next shape before we return
                f.seek(next)
                return None
        # Shape types with parts
        if shapeType in (3,5,13,15,23,25,31):
            nParts = unpack("<i", f.read(4))[0]
        # Shape types with points
        if shapeType in (3,5,8,13,15,18,23,25,28,31):
            nPoints = unpack("<i", f.read(4))[0]
        # Read parts
        if nParts:
            record.parts = _Array('i', unpack("<%si" % nParts, f.read(nParts * 4)))
        # Read part types for Multipatch - 31
        if shapeType == 31:
            record.partTypes = _Array('i', unpack("<%si" % nParts, f.read(nParts * 4)))
        # Read points - produces a list of [x,y] values
        if nPoints:
            flat = unpack("<%sd" % (2 * nPoints), f.read(16*nPoints))
            record.points = list(izip(*(iter(flat),) * 2))
        # Read z extremes and values
        if shapeType in (13,15,18,31):
            (zmin, zmax) = unpack("<2d", f.read(16))
            record.z = _Array('d', unpack("<%sd" % nPoints, f.read(nPoints * 8)))
        # Read m extremes and values
        if shapeType in (13,15,18,23,25,28,31):
            if next - f.tell() >= 16:
                (mmin, mmax) = unpack("<2d", f.read(16))
            # Measure values less than -10e38 are nodata values according to the spec
            if next - f.tell() >= nPoints * 8:
                record.m = []
                for m in _Array('d', unpack("<%sd" % nPoints, f.read(nPoints * 8))):
                    if m > NODATA:
                        record.m.append(m)
                    else:
                        record.m.append(None)
            else:
                record.m = [None for _ in range(nPoints)]
        # Read a single point
        if shapeType in (1,11,21):
            record.points = [_Array('d', unpack("<2d", f.read(16)))]
        # Read a single Z value
        if shapeType == 11:
            record.z = list(unpack("<d", f.read(8)))
        # Read a single M value
        if shapeType in (21,11):
            if next - f.tell() >= 8:
                (m,) = unpack("<d", f.read(8))
            else:
                m = NODATA
            # Measure values less than -10e38 are nodata values according to the spec
            if m > NODATA:
                record.m = [m]
            else:
                record.m = [None]
        # Seek to the end of this record as defined by the record header because
        # the shapefile spec doesn't require the actual content to meet the header
        # definition.  Probably allowed for lazy feature deletion. 
        f.seek(next)
        return record
    def iterRecords(self, fields=None):
        """Returns a generator of records in a dbf file.
        Useful for large shapefiles or dbf files.
        To only read some of the fields, specify the 'fields' arg as a
        list of one or more fieldnames.
        """
        if self.numRecords is None:
            self.__dbfHeader()
        f = self.__getFileObj(self.dbf)
        f.seek(self.__dbfHdrLength)
        fieldTuples,recLookup,recStruct = self.__recordFields(fields)
        for i in range(self.numRecords):
            r = self.__record(oid=i, fieldTuples=fieldTuples, recLookup=recLookup, recStruct=recStruct)
            if r:
                yield r
    @property 
    def __json__(self):
        shapeRecords=self.shapeRecords
        fcoll=shapeRecords.__json__
        fcoll['bbox']=list(self.bbox)
        return fcoll
    def __record(self, fieldTuples, recLookup, recStruct, oid=None):
        """Reads and returns a dbf record row as a list of values. Requires specifying
        a list of field info tuples 'fieldTuples', a record name-index dict 'recLookup', 
        and a Struct instance 'recStruct' for unpacking these fields. 
        """
        f = self.__getFileObj(self.dbf)

        recordContents = recStruct.unpack(f.read(recStruct.size))
        
        # deletion flag field is always unpacked as first value (see __recordFmt)
        if recordContents[0] != b' ':
            # deleted record
            return None

        # drop deletion flag from values
        recordContents = recordContents[1:]

        # check that values match fields
        if len(fieldTuples) != len(recordContents):
            raise Exception('Number of record values ({}) is different from the requested \
                            number of fields ({})'.format(len(recordContents), len(fieldTuples)))

        # parse each value
        record = []
        for (name, typ, size, deci),value in zip(fieldTuples, recordContents):
            if typ in ("N","F"):
                # numeric or float: number stored as a string, right justified, and padded with blanks to the width of the field. 
                value = value.split(b'\0')[0]
                value = value.replace(b'*', b'')  # QGIS NULL is all '*' chars
                if value == b'':
                    value = None
                elif deci:
                    try:
                        value = float(value)
                    except ValueError:
                        #not parseable as float, set to None
                        value = None
                else:
                    # force to int
                    try:
                        # first try to force directly to int.
                        # forcing a large int to float and back to int
                        # will lose information and result in wrong nr.
                        value = int(value) 
                    except ValueError:
                        # forcing directly to int failed, so was probably a float.
                        try:
                            value = int(float(value))
                        except ValueError:
                            #not parseable as int, set to None
                            value = None
            elif typ == 'D':
                # date: 8 bytes - date stored as a string in the format YYYYMMDD.
                if not value.replace(b'\x00', b'').replace(b' ', b'').replace(b'0', b''):
                    # dbf date field has no official null value
                    # but can check for all hex null-chars, all spaces, or all 0s (QGIS null)
                    value = None
                else:
                    try:
                        # return as python date object
                        y, m, d = int(value[:4]), int(value[4:6]), int(value[6:8])
                        from datetime import date
                        value = date(y, m, d)
                    except:
                        # if invalid date, just return as unicode string so user can decide
                        value = u(value.strip())
            elif typ == 'L':
                # logical: 1 byte - initialized to 0x20 (space) otherwise T or F.
                if value == b" ":
                    value = None # space means missing or not yet set
                else:
                    if value in b'YyTt1':
                        value = True
                    elif value in b'NnFf0':
                        value = False
                    else:
                        value = None # unknown value is set to missing
            else:
                # anything else is forced to string/unicode
                value = u(value, self.encoding, self.encodingErrors)
                value = value.strip().rstrip('\x00') # remove null-padding at end of strings
            record.append(value)

        return _Record(recLookup, record, oid)



    

class ShapeRecords(list):
    @property
    def __json__(self):
        collection={'type':'FeatureCollection','features':[shpR.__json__ for shpR in self]}
        return collection;




import array
class _Array(array.array):
    def __repr__(self):
        return str(self.tolist());
    




class Shape(object):
    def __init__(self, shapeType=NULL, points=None, parts=None, partTypes=None, oid=None):
        """Stores the geometry of the different shape types
        specified in the Shapefile spec. Shape types are
        usually point, polyline, or polygons. Every shape type
        except the "Null" type contains points at some level for
        example vertices in a polygon. If a shape type has
        multiple shapes containing points within a single
        geometry record then those shapes are called parts. Parts
        are designated by their starting index in geometry record's
        list of shapes. For MultiPatch geometry, partTypes designates
        the patch type of each of the parts. 
        """
        self.shapeType = shapeType
        self.points = points or []
        self.parts = parts or []
        if partTypes:
            self.partTypes = partTypes
        
        # and a dict to silently record any errors encountered
        self._errors = {}
        
        # add oid
        if oid is not None:
            self.__oid = oid
        else:
            self.__oid = -1

    @property
    def __geo_interface__(self):
        if self.shapeType in [POINT, POINTM, POINTZ]:
            # point
            if len(self.points) == 0:
                # the shape has no coordinate information, i.e. is 'empty'
                # the geojson spec does not define a proper null-geometry type
                # however, it does allow geometry types with 'empty' coordinates to be interpreted as null-geometries
                return {'type':'Point', 'coordinates':tuple()}
            else:
                return {
                'type': 'Point',
                'coordinates': tuple(self.points[0])
                }
        elif self.shapeType in [MULTIPOINT, MULTIPOINTM, MULTIPOINTZ]:
            if len(self.points) == 0:
                # the shape has no coordinate information, i.e. is 'empty'
                # the geojson spec does not define a proper null-geometry type
                # however, it does allow geometry types with 'empty' coordinates to be interpreted as null-geometries
                return {'type':'MultiPoint', 'coordinates':[]}
            else:
                # multipoint
                return {
                'type': 'MultiPoint',
                'coordinates': [tuple(p) for p in self.points]
                }
        elif self.shapeType in [POLYLINE, POLYLINEM, POLYLINEZ]:
            if len(self.parts) == 0:
                # the shape has no coordinate information, i.e. is 'empty'
                # the geojson spec does not define a proper null-geometry type
                # however, it does allow geometry types with 'empty' coordinates to be interpreted as null-geometries
                return {'type':'LineString', 'coordinates':[]}
            elif len(self.parts) == 1:
                # linestring
                return {
                'type': 'LineString',
                'coordinates': [tuple(p) for p in self.points]
                }
            else:
                # multilinestring
                ps = None
                coordinates = []
                for part in self.parts:
                    if ps == None:
                        ps = part
                        continue
                    else:
                        coordinates.append([tuple(p) for p in self.points[ps:part]])
                        ps = part
                else:
                    coordinates.append([tuple(p) for p in self.points[part:]])
                return {
                'type': 'MultiLineString',
                'coordinates': coordinates
                }
        elif self.shapeType in [POLYGON, POLYGONM, POLYGONZ]:
            if len(self.parts) == 0:
                # the shape has no coordinate information, i.e. is 'empty'
                # the geojson spec does not define a proper null-geometry type
                # however, it does allow geometry types with 'empty' coordinates to be interpreted as null-geometries
                return {'type':'Polygon', 'coordinates':[]}
            else:
                # get all polygon rings
                rings = []
                for i in xrange(len(self.parts)):
                    # get indexes of start and end points of the ring
                    start = self.parts[i]
                    try:
                        end = self.parts[i+1]
                    except IndexError:
                        end = len(self.points)

                    # extract the points that make up the ring
                    ring = [tuple(p) for p in self.points[start:end]]
                    rings.append(ring)

                # organize rings into list of polygons, where each polygon is defined as list of rings.
                # the first ring is the exterior and any remaining rings are holes (same as GeoJSON). 
                polys = organize_polygon_rings(rings, self._errors)
                
                # if VERBOSE is True, issue detailed warning about any shape errors
                # encountered during the Shapefile to GeoJSON conversion
                VERBOSE=False;
                if VERBOSE and self._errors: 
                    header = 'Possible issue encountered when converting Shape #{} to GeoJSON: '.format(self.oid)
                    orphans = self._errors.get('polygon_orphaned_holes', None)
                    if orphans:
                        msg = header + 'Shapefile format requires that all polygon interior holes be contained by an exterior ring, \
but the Shape contained interior holes (defined by counter-clockwise orientation in the shapefile format) that were \
orphaned, i.e. not contained by any exterior rings. The rings were still included but were \
encoded as GeoJSON exterior rings instead of holes.'
                        logger.warning(msg)
                    only_holes = self._errors.get('polygon_only_holes', None)
                    if only_holes:
                        msg = header + 'Shapefile format requires that polygons contain at least one exterior ring, \
but the Shape was entirely made up of interior holes (defined by counter-clockwise orientation in the shapefile format). The rings were \
still included but were encoded as GeoJSON exterior rings instead of holes.'
                        logger.warning(msg)

                # return as geojson
                if len(polys) == 1:
                    return {
                    'type': 'Polygon',
                    'coordinates': polys[0]
                    }
                else:
                    return {
                    'type': 'MultiPolygon',
                    'coordinates': polys
                    }

        else:
            raise Exception('Shape type "%s" cannot be represented as GeoJSON.' % SHAPETYPE_LOOKUP[self.shapeType])

    @staticmethod
    def _from_geojson(geoj):
        # create empty shape
        shape = Shape()
        # set shapeType
        geojType = geoj["type"] if geoj else "Null"
        if geojType == "Null":
            shapeType = NULL
        elif geojType == "Point":
            shapeType = POINT
        elif geojType == "LineString":
            shapeType = POLYLINE
        elif geojType == "Polygon":
            shapeType = POLYGON
        elif geojType == "MultiPoint":
            shapeType = MULTIPOINT
        elif geojType == "MultiLineString":
            shapeType = POLYLINE
        elif geojType == "MultiPolygon":
            shapeType = POLYGON
        else:
            raise Exception("Cannot create Shape from GeoJSON type '%s'" % geojType)
        shape.shapeType = shapeType
        
        # set points and parts
        if geojType == "Point":
            shape.points = [ geoj["coordinates"] ]
            shape.parts = [0]
        elif geojType in ("MultiPoint","LineString"):
            shape.points = geoj["coordinates"]
            shape.parts = [0]
        elif geojType in ("Polygon"):
            points = []
            parts = []
            index = 0
            for i,ext_or_hole in enumerate(geoj["coordinates"]):
                # although the latest GeoJSON spec states that exterior rings should have 
                # counter-clockwise orientation, we explicitly check orientation since older 
                # GeoJSONs might not enforce this. 
                if i == 0 and not is_cw(ext_or_hole):
                    # flip exterior direction
                    ext_or_hole = rewind(ext_or_hole)
                elif i > 0 and is_cw(ext_or_hole):
                    # flip hole direction
                    ext_or_hole = rewind(ext_or_hole)
                points.extend(ext_or_hole)
                parts.append(index)
                index += len(ext_or_hole)
            shape.points = points
            shape.parts = parts
        elif geojType in ("MultiLineString"):
            points = []
            parts = []
            index = 0
            for linestring in geoj["coordinates"]:
                points.extend(linestring)
                parts.append(index)
                index += len(linestring)
            shape.points = points
            shape.parts = parts
        elif geojType in ("MultiPolygon"):
            points = []
            parts = []
            index = 0
            for polygon in geoj["coordinates"]:
                for i,ext_or_hole in enumerate(polygon):
                    # although the latest GeoJSON spec states that exterior rings should have 
                    # counter-clockwise orientation, we explicitly check orientation since older 
                    # GeoJSONs might not enforce this. 
                    if i == 0 and not is_cw(ext_or_hole):
                        # flip exterior direction
                        ext_or_hole = rewind(ext_or_hole)
                    elif i > 0 and is_cw(ext_or_hole):
                        # flip hole direction
                        ext_or_hole = rewind(ext_or_hole)
                    points.extend(ext_or_hole)
                    parts.append(index)
                    index += len(ext_or_hole)
            shape.points = points
            shape.parts = parts
        return shape

    @property
    def oid(self):
        """The index position of the shape in the original shapefile"""
        return self.__oid

    @property
    def shapeTypeName(self):
        return SHAPETYPE_LOOKUP[self.shapeType]

    def __repr__(self):
        return 'Shape #{}: {}'.format(self.__oid, self.shapeTypeName)
    

class _Record(list):
    """
    A class to hold a record. Subclasses list to ensure compatibility with
    former work and to reuse all the optimizations of the builtin list.
    In addition to the list interface, the values of the record
    can also be retrieved using the field's name. For example if the dbf contains
    a field ID at position 0, the ID can be retrieved with the position, the field name
    as a key, or the field name as an attribute.

    >>> # Create a Record with one field, normally the record is created by the Reader class
    >>> r = _Record({'ID': 0}, [0])
    >>> print(r[0])
    >>> print(r['ID'])
    >>> print(r.ID)
    """

    def __init__(self, field_positions, values, oid=None):
        """
        A Record should be created by the Reader class

        :param field_positions: A dict mapping field names to field positions
        :param values: A sequence of values
        :param oid: The object id, an int (optional)
        """
        self.__field_positions = field_positions
        if oid is not None:
            self.__oid = oid
        else:
            self.__oid = -1
        list.__init__(self, values)

    def __getattr__(self, item):
        """
        __getattr__ is called if an attribute is used that does
        not exist in the normal sense. For example r=Record(...), r.ID
        calls r.__getattr__('ID'), but r.index(5) calls list.index(r, 5)
        :param item: The field name, used as attribute
        :return: Value of the field
        :raises: AttributeError, if item is not a field of the shapefile
                and IndexError, if the field exists but the field's 
                corresponding value in the Record does not exist
        """
        try:
            index = self.__field_positions[item]
            return list.__getitem__(self, index)
        except KeyError:
            raise AttributeError('{} is not a field name'.format(item))
        except IndexError:
            raise IndexError('{} found as a field but not enough values available.'.format(item))

    def __setattr__(self, key, value):
        """
        Sets a value of a field attribute
        :param key: The field name
        :param value: the value of that field
        :return: None
        :raises: AttributeError, if key is not a field of the shapefile
        """
        if key.startswith('_'):  # Prevent infinite loop when setting mangled attribute
            return list.__setattr__(self, key, value)
        try:
            index = self.__field_positions[key]
            return list.__setitem__(self, index, value)
        except KeyError:
            raise AttributeError('{} is not a field name'.format(key))

    def __getitem__(self, item):
        """
        Extends the normal list item access with
        access using a fieldname

        For example r['ID'], r[0]
        :param item: Either the position of the value or the name of a field
        :return: the value of the field
        """
        try:
            return list.__getitem__(self, item)
        except TypeError:
            try:
                index = self.__field_positions[item]
            except KeyError:
                index = None
        if index is not None:
            return list.__getitem__(self, index)
        else:
            raise IndexError('"{}" is not a field name and not an int'.format(item))

    def __setitem__(self, key, value):
        """
        Extends the normal list item access with
        access using a fieldname

        For example r['ID']=2, r[0]=2
        :param key: Either the position of the value or the name of a field
        :param value: the new value of the field
        """
        try:
            return list.__setitem__(self, key, value)
        except TypeError:
            index = self.__field_positions.get(key)
            if index is not None:
                return list.__setitem__(self, index, value)
            else:
                raise IndexError('{} is not a field name and not an int'.format(key))

    @property
    def oid(self):
        """The index position of the record in the original shapefile"""
        return self.__oid

    def as_dict(self, date_strings=False):
        """
        Returns this Record as a dictionary using the field names as keys
        :return: dict
        """
        from datetime import date
        dct = dict((f, self[i]) for f, i in self.__field_positions.items())
        if date_strings:
            for k,v in dct.items():
                if isinstance(v, date):
                    dct[k] = '{:04d}{:02d}{:02d}'.format(v.year, v.month, v.day)
        return dct

    def __repr__(self):
        return 'Record #{}: {}'.format(self.__oid, list(self))

    def __dir__(self):
        """
        Helps to show the field names in an interactive environment like IPython.
        See: http://ipython.readthedocs.io/en/stable/config/integrating.html

        :return: List of method names and fields
        """
        default = list(dir(type(self))) # default list methods and attributes of this class
        fnames = list(self.__field_positions.keys()) # plus field names (random order if Python version < 3.6)
        return default + fnames

class ShapeRecord(object):
    """A ShapeRecord object containing a shape along with its attributes.
    Provides the GeoJSON __geo_interface__ to return a Feature dictionary."""
    def __init__(self, shape=None, record=None):
        self.shape = shape
        self.record = record

    @property
    def __json__(self):
        return {'type': 'Feature',
                'properties': self.record.as_dict(date_strings=True),
                'geometry': None if self.shape.shapeType == NULL else self.shape.__geo_interface__}

def organize_polygon_rings(rings, return_errors=None):
    '''Organize a list of coordinate rings into one or more polygons with holes.
    Returns a list of polygons, where each polygon is composed of a single exterior
    ring, and one or more interior holes. If a return_errors dict is provided (optional), 
    any errors encountered will be added to it. 

    Rings must be closed, and cannot intersect each other (non-self-intersecting polygon).
    Rings are determined as exteriors if they run in clockwise direction, or interior
    holes if they run in counter-clockwise direction. This method is used to construct
    GeoJSON (multi)polygons from the shapefile polygon shape type, which does not
    explicitly store the structure of the polygons beyond exterior/interior ring orientation. 
    '''
    # first iterate rings and classify as exterior or hole
    exteriors = []
    holes = []
    for ring in rings:
        # shapefile format defines a polygon as a sequence of rings
        # where exterior rings are clockwise, and holes counterclockwise
        if is_cw(ring):
            # ring is exterior
            exteriors.append(ring)
        else:
            # ring is a hole
            holes.append(ring)
                
    # if only one exterior, then all holes belong to that exterior
    if len(exteriors) == 1:
        # exit early
        poly = [exteriors[0]] + holes
        polys = [poly]
        return polys

    # multiple exteriors, ie multi-polygon, have to group holes with correct exterior
    # shapefile format does not specify which holes belong to which exteriors
    # so have to do efficient multi-stage checking of hole-to-exterior containment
    elif len(exteriors) > 1:
        # exit early if no holes
        if not holes:
            polys = []
            for ext in exteriors:
                poly = [ext]
                polys.append(poly)
            return polys
        
        # first determine each hole's candidate exteriors based on simple bbox contains test
        hole_exteriors = dict([(hole_i,[]) for hole_i in xrange(len(holes))])
        exterior_bboxes = [ring_bbox(ring) for ring in exteriors]
        for hole_i in hole_exteriors.keys():
            hole_bbox = ring_bbox(holes[hole_i])
            for ext_i,ext_bbox in enumerate(exterior_bboxes):
                if bbox_contains(ext_bbox, hole_bbox):
                    hole_exteriors[hole_i].append( ext_i )

        # then, for holes with still more than one possible exterior, do more detailed hole-in-ring test
        for hole_i,exterior_candidates in hole_exteriors.items():
            
            if len(exterior_candidates) > 1:
                # get hole sample point
                ccw = not is_cw(holes[hole_i])
                hole_sample = ring_sample(holes[hole_i], ccw=ccw)
                # collect new exterior candidates
                new_exterior_candidates = []
                for ext_i in exterior_candidates:
                    # check that hole sample point is inside exterior
                    hole_in_exterior = ring_contains_point(exteriors[ext_i], hole_sample)
                    if hole_in_exterior:
                        new_exterior_candidates.append(ext_i)

                # set new exterior candidates
                hole_exteriors[hole_i] = new_exterior_candidates

        # if still holes with more than one possible exterior, means we have an exterior hole nested inside another exterior's hole
        for hole_i,exterior_candidates in hole_exteriors.items():
            
            if len(exterior_candidates) > 1:
                # exterior candidate with the smallest area is the hole's most immediate parent
                ext_i = sorted(exterior_candidates, key=lambda x: abs(signed_area(exteriors[x], fast=True)))[0]
                hole_exteriors[hole_i] = [ext_i]

        # separate out holes that are orphaned (not contained by any exterior)
        orphan_holes = []
        for hole_i,exterior_candidates in list(hole_exteriors.items()):
            if not exterior_candidates:
                orphan_holes.append( hole_i )
                del hole_exteriors[hole_i]
                continue

        # each hole should now only belong to one exterior, group into exterior-holes polygons
        polys = []
        for ext_i,ext in enumerate(exteriors):
            poly = [ext]
            # find relevant holes
            poly_holes = []
            for hole_i,exterior_candidates in list(hole_exteriors.items()):
                # hole is relevant if previously matched with this exterior
                if exterior_candidates[0] == ext_i:
                    poly_holes.append( holes[hole_i] )
            poly += poly_holes
            polys.append(poly)

        # add orphan holes as exteriors
        for hole_i in orphan_holes:
            ext = holes[hole_i]
            # add as single exterior without any holes
            poly = [ext]
            polys.append(poly)

        if orphan_holes and return_errors is not None:
            return_errors['polygon_orphaned_holes'] = len(orphan_holes)

        return polys

    # no exteriors, be nice and assume due to incorrect winding order
    else:
        if return_errors is not None:
            return_errors['polygon_only_holes'] = len(holes)
        exteriors = holes
        # add as single exterior without any holes
        polys = [[ext] for ext in exteriors]
        return polys
    




def is_cw(coords):
    """Returns True if a polygon ring has clockwise orientation, determined
    by a negatively signed area. 
    """
    area2 = signed_area(coords, fast=True)
    return area2 < 0





def signed_area(coords, fast=False):
    """Return the signed area enclosed by a ring using the linear time
    algorithm. A value >= 0 indicates a counter-clockwise oriented ring.
    A faster version is possible by setting 'fast' to True, which returns
    2x the area, e.g. if you're only interested in the sign of the area.
    """
    xs, ys = map(list, list(zip(*coords))[:2]) # ignore any z or m values
    xs.append(xs[1])
    ys.append(ys[1])
    area2 = sum(xs[i]*(ys[i+1]-ys[i-1]) for i in range(1, len(coords)))
    if fast:
        return area2
    else:
        return area2 / 2.0
    


def rewind(coords):
    """Returns the input coords in reversed order.
    """
    return list(reversed(coords))

def ring_bbox(coords):
    """Calculates and returns the bounding box of a ring.
    """
    xs,ys = zip(*coords)
    bbox = min(xs),min(ys),max(xs),max(ys)
    return bbox

def bbox_overlap(bbox1, bbox2):
    """Tests whether two bounding boxes overlap, returning a boolean
    """
    xmin1,ymin1,xmax1,ymax1 = bbox1
    xmin2,ymin2,xmax2,ymax2 = bbox2
    overlap = (xmin1 <= xmax2 and xmax1 >= xmin2 and ymin1 <= ymax2 and ymax1 >= ymin2)
    return overlap

def bbox_contains(bbox1, bbox2):
    """Tests whether bbox1 fully contains bbox2, returning a boolean
    """
    xmin1,ymin1,xmax1,ymax1 = bbox1
    xmin2,ymin2,xmax2,ymax2 = bbox2
    contains = (xmin1 < xmin2 and xmax1 > xmax2 and ymin1 < ymin2 and ymax1 > ymax2)
    return contains

def ring_contains_point(coords, p):
    """Fast point-in-polygon crossings algorithm, MacMartin optimization.

    Adapted from code by Eric Haynes
    http://www.realtimerendering.com/resources/GraphicsGems//gemsiv/ptpoly_haines/ptinpoly.c
    
    Original description:
        Shoot a test ray along +X axis.  The strategy, from MacMartin, is to
        compare vertex Y values to the testing point's Y and quickly discard
        edges which are entirely to one side of the test ray.
    """
    tx,ty = p

    # get initial test bit for above/below X axis
    vtx0 = coords[0]
    yflag0 = ( vtx0[1] >= ty )

    inside_flag = False
    for vtx1 in coords[1:]: 
        yflag1 = ( vtx1[1] >= ty )
        # check if endpoints straddle (are on opposite sides) of X axis
        # (i.e. the Y's differ); if so, +X ray could intersect this edge.
        if yflag0 != yflag1: 
            xflag0 = ( vtx0[0] >= tx )
            # check if endpoints are on same side of the Y axis (i.e. X's
            # are the same); if so, it's easy to test if edge hits or misses.
            if xflag0 == ( vtx1[0] >= tx ):
                # if edge's X values both right of the point, must hit
                if xflag0:
                    inside_flag = not inside_flag
            else:
                # compute intersection of pgon segment with +X ray, note
                # if >= point's X; if so, the ray hits it.
                if ( vtx1[0] - (vtx1[1]-ty) * ( vtx0[0]-vtx1[0]) / (vtx0[1]-vtx1[1]) ) >= tx:
                    inside_flag = not inside_flag

        # move to next pair of vertices, retaining info as possible
        yflag0 = yflag1
        vtx0 = vtx1

    return inside_flag


def ring_sample(coords, ccw=False):
    """Return a sample point guaranteed to be within a ring, by efficiently
    finding the first centroid of a coordinate triplet whose orientation
    matches the orientation of the ring and passes the point-in-ring test.
    The orientation of the ring is assumed to be clockwise, unless ccw
    (counter-clockwise) is set to True. 
    """
    triplet = []
    def itercoords():
        # iterate full closed ring
        for p in coords:
            yield p
        # finally, yield the second coordinate to the end to allow checking the last triplet
        yield coords[1]
        
    for p in itercoords(): 
        # add point to triplet (but not if duplicate)
        if p not in triplet:
            triplet.append(p)
            
        # new triplet, try to get sample
        if len(triplet) == 3:
            # check that triplet does not form a straight line (not a triangle)
            is_straight_line = (triplet[0][1] - triplet[1][1]) * (triplet[0][0] - triplet[2][0]) == (triplet[0][1] - triplet[2][1]) * (triplet[0][0] - triplet[1][0])
            if not is_straight_line:
                # get triplet orientation
                closed_triplet = triplet + [triplet[0]]
                triplet_ccw = not is_cw(closed_triplet)
                # check that triplet has the same orientation as the ring (means triangle is inside the ring)
                if ccw == triplet_ccw:
                    # get triplet centroid
                    xs,ys = zip(*triplet)
                    xmean,ymean = sum(xs) / 3.0, sum(ys) / 3.0
                    # check that triplet centroid is truly inside the ring
                    if ring_contains_point(coords, (xmean,ymean)):
                        return xmean,ymean

            # failed to get sample point from this triplet
            # remove oldest triplet coord to allow iterating to next triplet
            triplet.pop(0)
            
    else:
        raise Exception('Unexpected error: Unable to find a ring sample point.')