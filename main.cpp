#include <math.h>
#include <iostream>

// OpenCASCADE header files
#include <BRep_Tool.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <Poly_Triangulation.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <STEPControl_Reader.hxx>
#include <NCollection_List.hxx>
#include <HLRBRep_Algo.hxx>
#include <gp_Ax2.hxx>
#include <HLRBRep_PolyAlgo.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <HLRBRep_HLRToShape.hxx>
#include <STEPControl_Writer.hxx>
#include <BRep_TEdge.hxx>
#include <Geom2d_BSplineCurve.hxx>

////libdxf header files
//#include <spline.h>

// Prototypes
void shape_bounding_box(TopoDS_Shape &s, double &xmin, double &xmax, double &ymin, double &ymax, double &zmin, double &zmax);

/*
 * (c) allpccloud GmbH, 2022
 * Author: Markus Finck
 *
 * This test program reads a STEP file (3d CAD geometry), creates a triangle mesh
 * and computes the bounding box using the mesh nodes.
 *
 * The OpenCASCADE library is used as the 3d CAD kernel (version 6.9.0 or higher).
 *
 * Compiler GNU C++ (gcc 7.5.0).
 * Please see the file SIGTestProject.pro for compiler and linker flags.
 *
 */
int main(int argc, char **argv)
{
    STEPControl_Reader reader;

    // Check command line
    if(argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " step-file" << std::endl;
        return 0;
    }

    // Read STEP file
    std::cout << "Read STEP file " << argv[1] << "..." << std::endl;
    reader.ReadFile(argv[1]);
    reader.TransferRoots();
    int nb_shapes = reader.NbShapes();

    // Check if shapes have not been read
    if(nb_shapes == 0)
    {
        std::cout << "STEP file " << argv[1] << " did not contain shapes" << std::endl;
        return 0;
    }

    // Shape(s) found - read in
    TopoDS_Shape shape = reader.OneShape();

    // Compute the bounding box
    double xmin, xmax, ymin, ymax, zmin, zmax;
    shape_bounding_box(shape, xmin, xmax, ymin, ymax, zmin, zmax);

    std::cout << "Bounding box:" << std::endl;
    std::cout << "xmin/max = " << xmin << ", " << xmax << std::endl;
    std::cout << "ymin/max = " << ymin << ", " << ymax << std::endl;
    std::cout << "zmin/max = " << zmin << ", " << zmax << std::endl;


    /*
    Pasted code 
    */
   // HLRAlgo_Projector 
    HLRAlgo_Projector myProjector(gp_Ax2(gp_Pnt(xmin, ymin, zmin), gp_Dir(0, 0, 1)));

    HLRBRep_Algo *myAlgo = new HLRBRep_Algo();
    // Add Shapes into the algorithm 
    //for (TopTools_ListOfShape::Iterator anIt(shape); anIt.More(); anIt.Next())
    //{
    //    myPolyAlgo->Load(anIt.Value());
    //}
    myAlgo->Add(shape);
    myAlgo->Projector(myProjector);

    // Build HLR 
    myAlgo->Update();

    // Set The Edge Status 
    myAlgo->Hide();

    // Build the extraction object : 
    HLRBRep_HLRToShape aHLRToShape(myAlgo);

    // extract the results : 
    TopoDS_Shape VCompound = aHLRToShape.VCompound();
    TopoDS_Shape Rg1LineVCompound =
        aHLRToShape.Rg1LineVCompound();
    TopoDS_Shape RgNLineVCompound =
        aHLRToShape.RgNLineVCompound();
    TopoDS_Shape OutLineVCompound =
        aHLRToShape.OutLineVCompound();
    TopoDS_Shape IsoLineVCompound =
        aHLRToShape.IsoLineVCompound();
    TopoDS_Shape HCompound = aHLRToShape.HCompound();
    TopoDS_Shape Rg1LineHCompound =
        aHLRToShape.Rg1LineHCompound();
    TopoDS_Shape RgNLineHCompound =
        aHLRToShape.RgNLineHCompound();
    TopoDS_Shape OutLineHCompound =
        aHLRToShape.OutLineHCompound();
    TopoDS_Shape IsoLineHCompound =
        aHLRToShape.IsoLineHCompound();

    if (VCompound.ShapeType() == TopAbs_ShapeEnum::TopAbs_COMPOUND)
    {
        std::cout << "VCompound.NbChildren:" << VCompound.NbChildren() <<std::endl;
        ////--------------------------------------
        //DxfFile fil = dxf_file_struct();
        //char s[14];
        //strcpy(s, "Test_file.dxf");
        //fil.fp = fopen(s, "w+");
        //fil.filename = s;
        ////--------------------------------------


        for (TopExp_Explorer exp(VCompound, TopAbs_ShapeEnum::TopAbs_EDGE);exp.More();exp.Next())
        {
            TopoDS_Edge aCurrentEdge = TopoDS::Edge(exp.Current());
            // do whatever with this edge
            Handle(BRep_TEdge) brep_TEdge = Handle(BRep_TEdge)::DownCast(aCurrentEdge.TShape());

            if (brep_TEdge.IsNull() == false)
            {
                BRep_ListOfCurveRepresentation curvesList = brep_TEdge->Curves();
                for (BRep_ListIteratorOfListOfCurveRepresentation iter(curvesList);iter.More();iter.Next())
                {
                    Handle(BRep_CurveRepresentation) curveRepr = iter.Value();
                    
                    
                    std::cout << "IsCurveOnSurface:" << curveRepr->IsCurveOnSurface() << std::endl;
                    std::cout << "IsPolygonOnClosedSurface:" << curveRepr->IsPolygonOnClosedSurface() << std::endl;
                    std::cout << "IsPolygonOnSurface:" << curveRepr->IsPolygonOnSurface() << std::endl;
                    std::cout << "IsPolygonOnClosedTriangulation:" << curveRepr->IsPolygonOnClosedTriangulation() << std::endl;
                    std::cout << "IsPolygonOnTriangulation:" << curveRepr->IsPolygonOnTriangulation() << std::endl;
                    std::cout << "IsCurve3D:" << curveRepr->IsCurve3D() << std::endl;
                    std::cout << "IsRegularity:" << curveRepr->IsRegularity() << std::endl;
                    std::cout << "IsCurveOnClosedSurface:" << curveRepr->IsCurveOnClosedSurface() << std::endl;
                    std::cout << "IsPolygon3D:" << curveRepr->IsPolygon3D() << std::endl;

                    if (curveRepr->IsCurveOnSurface())
                    {
                        Handle(Geom2d_Curve) curve = curveRepr->PCurve();
                        if (curve.IsNull() == false)
                        {
                            std::cout << "              curve" <<  std::endl;
                            std::cout << "              curve:" << curve->GetRefCount() <<std::endl;
                            std::cout << "              curve:" << curve->DynamicType() << std::endl;
                            std::cout << "              curve->DynamicType()->Name():"<< curve->DynamicType()->Name() << std::endl;
                            std::cout << "              curve:" << curve->get_type_descriptor() << std::endl;
                            std::cout << "              curve:" << curve->get_type_name() << std::endl;

                            TCollection_AsciiString str("Geom2d_BSplineCurve");
                            if (str.Search(curve->DynamicType()->Name())>-1)
                            {
                                Handle(Geom2d_BSplineCurve) bSplineCurve = Handle(Geom2d_BSplineCurve)::DownCast(curve);
                                std::cout << "             !!!!Geom2d_BSplineCurve" << std::endl;

                                Standard_Real UStart = bSplineCurve->FirstParameter();
                                Standard_Real ULast  = bSplineCurve->LastParameter();
                                std::cout << "             UStart:" << UStart << " ULast:" << ULast << std::endl;

                                Standard_Integer UKnotFirstIndex = bSplineCurve->FirstUKnotIndex();
                                Standard_Integer UKnotLastIndex = bSplineCurve->LastUKnotIndex();
                                std::cout << "             UKnotFirstIndex:" << UKnotFirstIndex << " UKnotLastIndex:" << UKnotLastIndex << std::endl;
                                
                                Standard_Integer NbKnots = bSplineCurve->NbKnots();
                                Standard_Integer NbPoles = bSplineCurve->NbPoles();
                                std::cout << "             NbKnots:" << NbKnots << " NbPoles:" << NbPoles << std::endl;
                                std::cout << std::endl;

                                gp_Pnt2d P0;
                                gp_Vec2d V;
                                Standard_Real U = UStart;
                                std::cout << "----------------------------" << std::endl;
                                bSplineCurve->D0(U, P0);
                                std::cout << "             D0 U:" << U << " PX:" << P0.X() << " PY:" << P0.Y() << std::endl;
          
                                bSplineCurve->D1(U, P0, V);
                                std::cout << "             D1 U:" << U << " PX:" << P0.X() << " PY:" << P0.Y() << std::endl;
                                std::cout << "                VX:" << V.X() << " VY:" << V.Y() << std::endl;
                                std::cout << "----------------------------" << std::endl;
                                
                                U = (ULast - UStart)/2 + UStart + 0.025;
                                bSplineCurve->D0(U, P0);
                                std::cout << "             D0 U:" << U << " PX:" << P0.X() << " PY:" << P0.Y() << std::endl;

                                bSplineCurve->D1(U, P0, V);
                                std::cout << "             D1 U:" << U << " PX:" << P0.X() << " PY:" << P0.Y() << std::endl;
                                std::cout << "                VX:" << V.X() << " VY:" << V.Y() << std::endl;
                                std::cout << "----------------------------" << std::endl;
                                
                                U = ULast;
                                bSplineCurve->D0(U, P0);
                                std::cout << "             D0 U:" << U << " PX:" << P0.X() << " PY:" << P0.Y() << std::endl;

                                bSplineCurve->D1(U, P0, V);
                                std::cout << "             D1 U:" << U << " PX:" << P0.X() << " PY:" << P0.Y() << std::endl;
                                std::cout << "                VX:" << V.X() << " VY:" << V.Y() << std::endl;
                                std::cout << "----------------------------" << std::endl;

                                //DxfSpline* spline = dxf_spline_new();

                                //DxfPoint* point = dxf_point_new();

                                //dxf_spline_write(&fil, spline);
                            }

                        }
                    }

                }
                std::cout << "----------------------" <<  std::endl;
            }
            
        }
        //dxf_file_write_eof(fil);
    }

    //STEPControl_Writer writerVCompound;
    //writerVCompound.Transfer(VCompound, STEPControl_GeometricCurveSet);
    //writerVCompound.Write("VCompound.stp");

    //STEPControl_Writer writerRg1LineVCompound;
    //writerRg1LineVCompound.Transfer(Rg1LineVCompound, STEPControl_GeometricCurveSet);
    //writerRg1LineVCompound.Write("Rg1LineVCompound.stp");

    //STEPControl_Writer writerRgNLineVCompound;
    //writerRgNLineVCompound.Transfer(RgNLineVCompound, STEPControl_GeometricCurveSet);
    //writerRgNLineVCompound.Write("RgNLineVCompound.stp");

    //STEPControl_Writer writerOutLineVCompound;
    //writerOutLineVCompound.Transfer(OutLineVCompound, STEPControl_GeometricCurveSet);
    // writerOutLineVCompound.Write("OutLineVCompound.stp");

    //STEPControl_Writer writerIsoLineVCompound;
    //writerIsoLineVCompound.Transfer(IsoLineVCompound, STEPControl_GeometricCurveSet);
    //writerIsoLineVCompound.Write("IsoLineVCompound.stp");


    //STEPControl_Writer writerHCompound;
    //writerHCompound.Transfer(HCompound, STEPControl_GeometricCurveSet);
    //writerHCompound.Write("HCompound.stp");

    //STEPControl_Writer writerRg1LineHCompound;
    //writerRg1LineHCompound.Transfer(Rg1LineHCompound, STEPControl_GeometricCurveSet);
    //writerRg1LineHCompound.Write("Rg1LineHCompound.stp");


    //STEPControl_Writer writerRgNLineHCompound;
    //writerRgNLineHCompound.Transfer(RgNLineHCompound, STEPControl_GeometricCurveSet);
    //writerRgNLineHCompound.Write("RgNLineHCompound.stp");

    //STEPControl_Writer writerOutLineHCompound;
    //writerOutLineHCompound.Transfer(OutLineHCompound, STEPControl_GeometricCurveSet);
    //writerOutLineHCompound.Write("OutLineHCompound.stp");

    //STEPControl_Writer writerIsoLineHCompound;
    //writerIsoLineHCompound.Transfer(IsoLineHCompound, STEPControl_GeometricCurveSet);
    //writerIsoLineHCompound.Write("IsoLineHCompound.stp");


    /*
     *
     * Please insert here source code for the test project.
     *
     * Links for the OpenCASCADE hidden line toolkits:
     * https://dev.opencascade.org/doc/overview/html/occt_contribution__documentation.html
     * and
     * https://dev.opencascade.org/doc/overview/html/occt_user_guides__modeling_algos.html#occt_modalg_10
     *
     * Choose a projection plane according to the bounding box and generate the feature edge projection.
     *
     * Export projection to DXF format which can be read by most CAD programs.
     * You can use the following libraries:
     * https://www.coin3d.org/dime/html/
     * https://github.com/bert/libdxf
     *
     * (more libraries capable of writing dxf are available)
     *
     */



    return 0;
}

void shape_bounding_box(TopoDS_Shape &s, double &xmin, double &xmax, double &ymin, double &ymax, double &zmin, double &zmax)
{
    TopExp_Explorer ex;
    double deflection = 0.001;
    double x,y,z;

    xmin = 1e20;
    xmax = -1e20;
    ymin = 1e20;
    ymax = -1e20;
    zmin = 1e20;
    zmax = -1e20;

    // Start meshing
    for(ex.Init(s,TopAbs_FACE);ex.More();ex.Next())
    {
        // Create mesh (triangulation), assignment to face is done automatically
        TopoDS_Face face = TopoDS::Face(ex.Current());

        if(!face.IsNull())
        {
            BRepMesh_IncrementalMesh im(face,deflection);

            // Get the initial triangulation
            TopLoc_Location L;
            Handle(Poly_Triangulation) poly_tri = BRep_Tool::Triangulation(face,L);
            if(poly_tri.IsNull()) return;
            gp_Trsf transf = L;

            // Get node list
            long num_nodes = poly_tri->NbNodes();
            Handle(TColgp_HArray1OfPnt) pnt_list = poly_tri->MapNodeArray();

            // Loop face mesh nodes
            for(long i=0;i<num_nodes;i++)
            { 
                // Extract node coordinates
                gp_Pnt pntFromArr = pnt_list->Value(i+1);
                gp_Pnt pnt = gp_Pnt(pntFromArr.Coord());
                pnt.Transform(transf);

                x = pnt.X();
                y = pnt.Y();
                z = pnt.Z();

                // Update bounding box
                xmin = fmin(x, xmin);
                xmax = fmax(x, xmax);
                ymin = fmin(y, ymin);
                ymax = fmax(y, ymax);
                zmin = fmin(z, zmin);
                zmax = fmax(z, zmax);
            }
        }
    }
}
