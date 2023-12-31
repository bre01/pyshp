from struct import *
import os
import array

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
izip = zip
xrange = range
SHAPETYPE_LOOKUP = {
    0: "NULL",
    1: "POINT",
    3: "POLYLINE",
    5: "POLYGON",
    8: "MULTIPOINT",
    11: "POINTZ",
    13: "POLYLINEZ",
    15: "POLYGONZ",
    18: "MULTIPOINTZ",
    21: "POINTM",
    23: "POLYLINEM",
    25: "POLYGONM",
    28: "MULTIPOINTM",
    31: "MULTIPATCH",
}

TRIANGLE_STRIP = 0
TRIANGLE_FAN = 1
OUTER_RING = 2
INNER_RING = 3
FIRST_RING = 4
RING = 5

PARTTYPE_LOOKUP = {
    0: "TRIANGLE_STRIP",
    1: "TRIANGLE_FAN",
    2: "OUTER_RING",
    3: "INNER_RING",
    4: "FIRST_RING",
    5: "RING",
}





class Shp(object):
    def __init__(self, *args, **kwargs):
        self.path = args[0]

        self.shp = None
        self.shx = None
        self.dbf = None
        self._file_to_close = []
        self.shapeName = "Not specified"
        self._offsets = []
        self.shpLength = None
        self.numRecords = None
        self.numShapes = None
        self.fields = []
        self.__dbfHdrLength = 0
        self.__fieldLookup = {}
        self.encoding = kwargs.pop("encoding", "utf-8")
        self.encodingErrors = kwargs.pop("encodingErrors", "strict")
        self.load(self.path)

    def load(self, shapefile=None):
        (shapeName, ext) = os.path.splitext(shapefile)
        self.shapeName = shapeName
        self.load_shp(shapeName)
        #self.load_shx(shapeName)
        if not (self.shp or self.dbf):
            raise Exception("Unable to open dbf or shp")
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
            raise Exception(
                "Shapefile Reader requires a shapefile or file-like object. (no shx file found"
            )
        # File length (16-bit word * 2 = bytes) - header length
        shx.seek(24)
        shxRecordLength = (unpack(">i", shx.read(4))[0] * 2) - 100
        self.numShapes = shxRecordLength // 8

    def load_shp(self, shapefile_name):
        shp_ext = "shp"
        self.shp = open("%s.%s" % (shapefile_name, shp_ext), "rb")
        self._file_to_close.append(self.shp)

    def load_shx(self, file_name):
        self.shx = open("%s.%s" % (file_name, "shx"), "rb")
        self._file_to_close.append(self.shx)


    def __shpHeader(self):
        shp = self.shp
        shp.seek(24)
        # file length big endian
        self.shpLength = unpack(">i", shp.read(4))[0] * 2
        # from word to byte

        # shape type little endian
        shp.seek(32)
        self.shapeType = unpack("<i", shp.read(4))[0]
        self.bbox = _Array("d", unpack("<4d", shp.read(32)))
        self.zbox = _Array("d", unpack("<2d", shp.read(16)))
        self.mbox = []
        for m in _Array("d", unpack("<2d", shp.read(16))):
            if m > NODATA:
                self.mbox.append(m)
            else:
                self.mbox.append(None)
                


    def iterShapes(self, bbox=None):
        """Returns a generator of shapes in a shapefile. Useful
        for handling large shapefiles.
        To only read shapes within a given spatial region, specify the 'bbox'
        arg as a list or tuple of xmin,ymin,xmax,ymax.
        """
        shp = self.shp
        # Found shapefiles which report incorrect
        # shp file length in the header. Can't trust
        # that so we seek to the end of the file
        # and figure it out.
        shp.seek(0, 2)
        shpLength = shp.tell()
        shp.seek(100)

        if self.numShapes:
            # Iterate exactly the number of shapes from shx header
            for i in range(self.numShapes):
                # MAYBE: check if more left of file or exit early?
                shape = self.setShapeIndex(oid=i)
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
                shape = self.setShapeIndex(oid=i)
                pos = shp.tell()
                if shape:
                    yield shape
                i += 1
            # Entire shp file consumed
            # Update the number of shapes and list of offsets
            assert i == len(offsets)
            self.numShapes = i
            self._offsets = offsets
    
    def __getFileObj(self,fileN):
        return fileN

    def setShapeIndex(self, oid=None):
        """Returns the header info and geometry for a single shape."""
        f = self.__getFileObj(self.shp)
        single_shape = Shape(oid=oid)
        nParts = nPoints = zmin = zmax = mmin = mmax = None
        (recNum, recLength) = unpack(">2i", f.read(8))
        # Determine the start of the next record
        next = f.tell() + (2 * recLength)
        shapeType = unpack("<i", f.read(4))[0]
        single_shape.shapeType = shapeType
        # For Null shapes create an empty points list for consistency
        if shapeType == 0:
            single_shape.points = []
        # All shape types capable of having a bounding box
        elif shapeType in (3, 5, 8, 13, 15, 18, 23, 25, 28, 31):
            single_shape.bbox = _Array("d", unpack("<4d", f.read(32)))
            # if bbox specified and no overlap, skip this shape
        # Shape types with parts
        if shapeType in (3, 5, 13, 15, 23, 25, 31):
            nParts = unpack("<i", f.read(4))[0]
        # Shape types with points

        if shapeType in (3, 5, 8, 13, 15, 18, 23, 25, 28, 31):
            nPoints = unpack("<i", f.read(4))[0]
        # Read parts
        if nParts:
            single_shape.parts = _Array("i", unpack("<%si" % nParts, f.read(nParts * 4)))


        # Read part types for Multipatch - 31
        """
        if shapeType == 31:
            single_shape.partTypes = _Array("i", unpack("<%si" % nParts, f.read(nParts * 4))) 
        """

        # Read points - produces a list of [x,y] values
        if nPoints:
            flat = unpack("<%sd" % (2 * nPoints), f.read(16 * nPoints))
            single_shape.points = list(izip(*(iter(flat),) * 2))
        # Read z extremes and values
        if shapeType in (13, 15, 18, 31):
            (zmin, zmax) = unpack("<2d", f.read(16))
            single_shape.z = _Array("d", unpack("<%sd" % nPoints, f.read(nPoints * 8)))
        # Read m extremes and values
        if shapeType in (13, 15, 18, 23, 25, 28, 31):
            if next - f.tell() >= 16:
                (mmin, mmax) = unpack("<2d", f.read(16))
            # Measure values less than -10e38 are nodata values according to the spec
            if next - f.tell() >= nPoints * 8:
                single_shape.m = []
                for m in _Array("d", unpack("<%sd" % nPoints, f.read(nPoints * 8))):
                    if m > NODATA:
                        single_shape.m.append(m)
                    else:
                        single_shape.m.append(None)
            else:
                single_shape.m = [None for _ in range(nPoints)]
        # Read a single point
        if shapeType in (1, 11, 21):
            single_shape.points = [_Array("d", unpack("<2d", f.read(16)))]
         
        # Seek to the end of this record as defined by the record header because
        # the shapefile spec doesn't require the actual content to meet the header
        # definition.  Probably allowed for lazy feature deletion.
        f.seek(next)
        return single_shape



    @property
    def shapeRecords(self):
        return ShapeRecords(self.iterShapeRecords())

    def iterShapeRecords(self):

        # the number of shape and record should be the same
        # in a typical(correct) shapefile 
        # so we combine dbf attribute and shape into the a shapeRecord
        for shape in self.iterShapes():
            yield ShapeRecord(shape=shape, record=None)

    @property
    def __json__(self):
        res= self.shapeRecords.__json__
        res["bbox"] = list(self.bbox)
        return res



class ShapeRecords(list):
    @property
    def __json__(self):
        collection = {
            "type": "FeatureCollection",
            "features": [shpRe.__json__ for shpRe in self],
        }
        return collection


class _Array(array.array):
    def __repr__(self):
        return str(self.tolist())


class Shape(object):
    def __init__(
        self, shapeType=NULL, points=None, parts=None, partTypes=None, oid=None
    ):
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
    def __json__(self):
        if self.shapeType in [POINT]:
            # point
            if len(self.points) == 0:
                # the shape has no coordinate information, i.e. is 'empty'
                # the geojson spec does not define a proper null-geometry type
                # however, it does allow geometry types with 'empty' coordinates to be interpreted as null-geometries
                return {"type": "Point", "coordinates": tuple()}
            else:
                return {"type": "Point", "coordinates": tuple(self.points[0])}

        elif self.shapeType in [POLYLINE]:
            if len(self.parts) == 0:
                # the shape has no coordinate information, i.e. is 'empty'
                # the geojson spec does not define a proper null-geometry type
                # however, it does allow geometry types with 'empty' coordinates to be interpreted as null-geometries
                return {"type": "LineString", "coordinates": []}
            elif len(self.parts) == 1:
                # linestring
                return {
                    "type": "LineString",
                    "coordinates": [tuple(p) for p in self.points],
                }
        elif self.shapeType in [POLYGON]:
        
                # the shape has no coordinate information, i.e. is 'empty'
                # the geojson spec does not define a proper null-geometry type

                # however, it does allow geometry types with 'empty' coordinates to be interpreted as null-geometries
            if len(self.parts)==0:
                return {"type": "Polygon", "coordinates": []}
            else:
                rings=[]
                for i in range(len(self.parts)):
                    start=self.parts[i]
                    try:
                        end=self.parts[i+1] 
                    except IndexError:
                        print("IndexError")
                        end=len(self.points)
                        
                    theRing=[p for p in self.points[start:end]]
                    rings.append(theRing)


                # get all polygon rings
                # if VERBOSE is True, issue detailed warning about any shape errors
                # encountered during the Shapefile to GeoJSON conversion
                

                # return as geojson
                return {"type": "Polygon", "coordinates": rings}



        else:
            raise Exception(
                "bad shape type"
            )


    def oid(self):
        """The index position of the shape in the original shapefile"""
        return self.__oid

    @property
    def shapeTypeName(self):
        return SHAPETYPE_LOOKUP[self.shapeType]

    def __repr__(self):
        return "Shape #{}: {}".format(self.__oid, self.shapeTypeName)


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
            raise AttributeError("{} is not a field name".format(item))
        except IndexError:
            raise IndexError(
                "{} found as a field but not enough values available.".format(item)
            )

    def __setattr__(self, key, value):
        """
        Sets a value of a field attribute
        :param key: The field name
        :param value: the value of that field
        :return: None
        :raises: AttributeError, if key is not a field of the shapefile
        """
        if key.startswith("_"):  # Prevent infinite loop when setting mangled attribute
            return list.__setattr__(self, key, value)
        try:
            index = self.__field_positions[key]
            return list.__setitem__(self, index, value)
        except KeyError:
            raise AttributeError("{} is not a field name".format(key))

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
                raise IndexError("{} is not a field name and not an int".format(key))

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
            for k, v in dct.items():
                if isinstance(v, date):
                    dct[k] = "{:04d}{:02d}{:02d}".format(v.year, v.month, v.day)
        return dct

    def __repr__(self):
        return "Record #{}: {}".format(self.__oid, list(self))

    def __dir__(self):
        """
        Helps to show the field names in an interactive environment like IPython.
        See: http://ipython.readthedocs.io/en/stable/config/integrating.html

        :return: List of method names and fields
        """
        default = list(
            dir(type(self))
        )  # default list methods and attributes of this class
        fnames = list(
            self.__field_positions.keys()
        )  # plus field names (random order if Python version < 3.6)
        return default + fnames


class ShapeRecord(object):
    """A ShapeRecord object containing a shape along with its attributes.
    Provides the __json__ property to return a Feature dictionary."""

    def __init__(self, shape=None, record=None):
        self.shape = shape
        self.record = record

    # line 493
    @property
    def __json__(self):

        # we don't have properties now, so return properties as null 
        return {
            "type": "Feature",
            "properties": None,
            "geometry": self.shape.__json__
        }


def organize_polygon_rings(rings, return_errors=None):
    """Organize a list of coordinate rings into one or more polygons with holes.
    Returns a list of polygons, where each polygon is composed of a single exterior
    ring, and one or more interior holes. If a return_errors dict is provided (optional),
    any errors encountered will be added to it.

    Rings must be closed, and cannot intersect each other (non-self-intersecting polygon).
    Rings are determined as exteriors if they run in clockwise direction, or interior
    holes if they run in counter-clockwise direction. This method is used to construct
    GeoJSON (multi)polygons from the shapefile polygon shape type, which does not
    explicitly store the structure of the polygons beyond exterior/interior ring orientation.
    """
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
        hole_exteriors = dict([(hole_i, []) for hole_i in xrange(len(holes))])
        exterior_bboxes = [ring_bbox(ring) for ring in exteriors]
        for hole_i in hole_exteriors.keys():
            hole_bbox = ring_bbox(holes[hole_i])
            for ext_i, ext_bbox in enumerate(exterior_bboxes):
                if bbox_contains(ext_bbox, hole_bbox):
                    hole_exteriors[hole_i].append(ext_i)

        # then, for holes with still more than one possible exterior, do more detailed hole-in-ring test
        for hole_i, exterior_candidates in hole_exteriors.items():
            if len(exterior_candidates) > 1:
                # get hole sample point
                ccw = not is_cw(holes[hole_i])
                hole_sample = ring_sample(holes[hole_i], ccw=ccw)
                # collect new exterior candidates
                new_exterior_candidates = []
                for ext_i in exterior_candidates:
                    # check that hole sample point is inside exterior
                    hole_in_exterior = ring_contains_point(
                        exteriors[ext_i], hole_sample
                    )
                    if hole_in_exterior:
                        new_exterior_candidates.append(ext_i)

                # set new exterior candidates
                hole_exteriors[hole_i] = new_exterior_candidates

        # if still holes with more than one possible exterior, means we have an exterior hole nested inside another exterior's hole
        for hole_i, exterior_candidates in hole_exteriors.items():
            if len(exterior_candidates) > 1:
                # exterior candidate with the smallest area is the hole's most immediate parent
                ext_i = sorted(
                    exterior_candidates,
                    key=lambda x: abs(signed_area(exteriors[x], fast=True)),
                )[0]
                hole_exteriors[hole_i] = [ext_i]

        # separate out holes that are orphaned (not contained by any exterior)
        orphan_holes = []
        for hole_i, exterior_candidates in list(hole_exteriors.items()):
            if not exterior_candidates:
                orphan_holes.append(hole_i)
                del hole_exteriors[hole_i]
                continue

        # each hole should now only belong to one exterior, group into exterior-holes polygons
        polys = []
        for ext_i, ext in enumerate(exteriors):
            poly = [ext]
            # find relevant holes
            poly_holes = []
            for hole_i, exterior_candidates in list(hole_exteriors.items()):
                # hole is relevant if previously matched with this exterior
                if exterior_candidates[0] == ext_i:
                    poly_holes.append(holes[hole_i])
            poly += poly_holes
            polys.append(poly)

        # add orphan holes as exteriors
        for hole_i in orphan_holes:
            ext = holes[hole_i]
            # add as single exterior without any holes
            poly = [ext]
            polys.append(poly)

        if orphan_holes and return_errors is not None:
            return_errors["polygon_orphaned_holes"] = len(orphan_holes)

        return polys

    # no exteriors, be nice and assume due to incorrect winding order
    else:
        if return_errors is not None:
            return_errors["polygon_only_holes"] = len(holes)
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
    xs, ys = map(list, list(zip(*coords))[:2])  # ignore any z or m values
    xs.append(xs[1])
    ys.append(ys[1])
    area2 = sum(xs[i] * (ys[i + 1] - ys[i - 1]) for i in range(1, len(coords)))
    if fast:
        return area2
    else:
        return area2 / 2.0


def rewind(coords):
    """Returns the input coords in reversed order."""
    return list(reversed(coords))


def ring_bbox(coords):
    """Calculates and returns the bounding box of a ring."""
    xs, ys = zip(*coords)
    bbox = min(xs), min(ys), max(xs), max(ys)
    return bbox


def bbox_overlap(bbox1, bbox2):
    """Tests whether two bounding boxes overlap, returning a boolean"""
    xmin1, ymin1, xmax1, ymax1 = bbox1
    xmin2, ymin2, xmax2, ymax2 = bbox2
    overlap = xmin1 <= xmax2 and xmax1 >= xmin2 and ymin1 <= ymax2 and ymax1 >= ymin2
    return overlap


def bbox_contains(bbox1, bbox2):
    """Tests whether bbox1 fully contains bbox2, returning a boolean"""
    xmin1, ymin1, xmax1, ymax1 = bbox1
    xmin2, ymin2, xmax2, ymax2 = bbox2
    contains = xmin1 < xmin2 and xmax1 > xmax2 and ymin1 < ymin2 and ymax1 > ymax2
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
    tx, ty = p

    # get initial test bit for above/below X axis
    vtx0 = coords[0]
    yflag0 = vtx0[1] >= ty

    inside_flag = False
    for vtx1 in coords[1:]:
        yflag1 = vtx1[1] >= ty
        # check if endpoints straddle (are on opposite sides) of X axis
        # (i.e. the Y's differ); if so, +X ray could intersect this edge.
        if yflag0 != yflag1:
            xflag0 = vtx0[0] >= tx
            # check if endpoints are on same side of the Y axis (i.e. X's
            # are the same); if so, it's easy to test if edge hits or misses.
            if xflag0 == (vtx1[0] >= tx):
                # if edge's X values both right of the point, must hit
                if xflag0:
                    inside_flag = not inside_flag
            else:
                # compute intersection of pgon segment with +X ray, note
                # if >= point's X; if so, the ray hits it.
                if (
                    vtx1[0] - (vtx1[1] - ty) * (vtx0[0] - vtx1[0]) / (vtx0[1] - vtx1[1])
                ) >= tx:
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
            is_straight_line = (triplet[0][1] - triplet[1][1]) * (
                triplet[0][0] - triplet[2][0]
            ) == (triplet[0][1] - triplet[2][1]) * (triplet[0][0] - triplet[1][0])
            if not is_straight_line:
                # get triplet orientation
                closed_triplet = triplet + [triplet[0]]
                triplet_ccw = not is_cw(closed_triplet)
                # check that triplet has the same orientation as the ring (means triangle is inside the ring)
                if ccw == triplet_ccw:
                    # get triplet centroid
                    xs, ys = zip(*triplet)
                    xmean, ymean = sum(xs) / 3.0, sum(ys) / 3.0
                    # check that triplet centroid is truly inside the ring
                    if ring_contains_point(coords, (xmean, ymean)):
                        return xmean, ymean

            # failed to get sample point from this triplet
            # remove oldest triplet coord to allow iterating to next triplet
            triplet.pop(0)

    else:
        raise Exception("Unexpected error: Unable to find a ring sample point.")
