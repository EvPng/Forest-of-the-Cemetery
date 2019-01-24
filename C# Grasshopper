using System;
using System.Collections;
using System.Collections.Generic;

using Rhino;
using Rhino.Geometry;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;



/// <summary>
/// This class will be instantiated on demand by the Script component.
/// </summary>
public class Script_Instance : GH_ScriptInstance
{
#region Utility functions
  /// <summary>Print a String to the [Out] Parameter of the Script component.</summary>
  /// <param name="text">String to print.</param>
  private void Print(string text) { /* Implementation hidden. */ }
  /// <summary>Print a formatted String to the [Out] Parameter of the Script component.</summary>
  /// <param name="format">String format.</param>
  /// <param name="args">Formatting parameters.</param>
  private void Print(string format, params object[] args) { /* Implementation hidden. */ }
  /// <summary>Print useful information about an object instance to the [Out] Parameter of the Script component. </summary>
  /// <param name="obj">Object instance to parse.</param>
  private void Reflect(object obj) { /* Implementation hidden. */ }
  /// <summary>Print the signatures of all the overloads of a specific method to the [Out] Parameter of the Script component. </summary>
  /// <param name="obj">Object instance to parse.</param>
  private void Reflect(object obj, string method_name) { /* Implementation hidden. */ }
#endregion

#region Members
  /// <summary>Gets the current Rhino document.</summary>
  private readonly RhinoDoc RhinoDocument;
  /// <summary>Gets the Grasshopper document that owns this script.</summary>
  private readonly GH_Document GrasshopperDocument;
  /// <summary>Gets the Grasshopper script component that owns this script.</summary>
  private readonly IGH_Component Component;
  /// <summary>
  /// Gets the current iteration count. The first call to RunScript() is associated with Iteration==0.
  /// Any subsequent call within the same solution will increment the Iteration count.
  /// </summary>
  private readonly int Iteration;
#endregion

  /// <summary>
  /// This procedure contains the user code. Input parameters are provided as regular arguments,
  /// Output parameters as ref arguments. You don't have to assign output parameters,
  /// they will have a default value.
  /// </summary>
  private void RunScript(List<double> s, List<Vector3d> v, int ry, int rx, double sizeX, double sizeY, Point3d focus, ref object A)
  {
    //creater an empty list to add the lines for the output
    List<Line> lines = new List<Line>();
    //compute steps of grid so that we can calculate the physical locations of points from the row/column indices
    double dx = sizeX / (double) (rx - 1.0);
    double dy = sizeY / (double) (ry - 1.0);
    //TODO GLOBAL1 [optional]: here you
    //Global Variables
    //General Density
    int scatterX = 2;
    int scatterY = 2;
    //The Range of Cemetery
    double rangeXY = 20.0;
    //Hilltop will not have plants
    double rangeZ = 4.0;
    //Control how many plants are randomly excluded
    double rangeRandom = 0.7;
    //Decide whether a land is smooth and gentle
    double angleZ = 0.2;
    //Scale Factor for elevation, used in generate random numbers
    double randomZ = 0.2;
    //Control the range of the very high pine trees surrounded the cemetery
    double rangeBelt = 0.2;
    //Control the grass levels
    double grassBoundary = 0.8;
    //The Step depth and width of path
    double step = 0.1;
    double stepL = 0.1;
    //General Random Generator
    int seed = 3;
    Random random = new Random(seed);
    //Useful Unit Vector and default plane
    Vector3d unitZ = new Vector3d(0, 0, 1);
    Plane oriPlane = new Plane(new Point3d(0, 0, 0), unitZ);
    //END GLOBAL1
    //start looping through all possible index pairs in the grid ry x rx
    double minDist = 10000.0;
    int minI = 0;
    int minJ = 0;
    int maxHeightI = 0;
    int maxHeightJ = 0;
    double maxHeight = 0.0;
    double max = 0.0;
    // Compute the highest elevation
    for (int n = 0; n < s.Count; n++)
    {
      if(s[n] > maxHeight){maxHeight = s[n];}
    }
    // Main Loop; Make components
    for (int j = 0; j < ry; ++j) {
      for(int i = 0; i < rx; ++i) {
        //compute the physical point location on the XY plane
        Point3d p = new Point3d(i * dx, j * dy, 0.0);
        //compute the serial index --> pickup values from the lists s and v
        int k = j * rx + i;
        // find distance and the closet point to focus
        double distanceXY = Math.Pow(p.X - focus.X, 2) + Math.Pow(p.Y - focus.Y, 2);
        // find the nearest point to visitor center
        if (distanceXY < minDist)
        {
          minDist = distanceXY;
          minI = i;
          minJ = j;
        }
        // store the location of the highest point in the landscape
        if (s[k] > max)
        {
          max = s[k];
          maxHeightI = i;
          maxHeightJ = j;
        }
        //TODO GLOBAL2 : here you need to decide whether to call the makeComponent method or not according to your criteria and logic
        double randomP = random.Next((int) (s[k] * randomZ), (int) (s[k] + 1));
        double ranCross = random.Next(k, k + 4);
        if ( i % scatterX == 0
          && j % scatterY == 0
          && distanceXY >= rangeXY
          && randomP >= rangeRandom
          && s[k] <= rangeZ)
        {
          //END GLOBAL2
          List<Line> component = makeComponent(p, s[k], v[k], focus, maxHeight); //call the makeCompoennt method from below to get hold of the local wireframe geometry

          //TODO GLOBAL3: here you can post-process the list of lines that the makeComponent returned before adding them to the master list
          // Scale
          Point3d anchor = new Point3d(p.X, p.Y, s[k]);
          Transform scale3D;
          if(distanceXY <= rangeXY * (1 + rangeBelt))
          {
            scale3D = Transform.Scale(anchor, 10.0 / (s[k] + 0.2));
          }
          else
          {
            scale3D = Transform.Scale(anchor, 1.0 * ( Math.Pow(Math.Abs(s[k] - 2.8), 2) + 0.5 + 0.05 * Math.Abs(distanceXY * 0.01 - 2)));
          }
          // Rotate
          Transform rotate2D = Transform.Rotation(Math.PI * Vector3d.VectorAngle(p - focus, unitZ), anchor);
          for (int m = 0, len = component.Count; m < len; m++)
          {
            Line temp = component[m];
            temp.Transform(scale3D);
            temp.Transform(rotate2D);
            component[m] = temp;
          }
          //END GLOBAL3
          lines.AddRange(component);
        }
        else if ( Vector3d.VectorAngle(v[k], unitZ) < angleZ)
        {
          if(distanceXY < rangeXY && ranCross < (k + 0.5))
          {
            //Randomly generate little crosses among grass
            List<Line> component = makeCross(p, s[k], focus);

            lines.AddRange(component);
          }
          else if(s[k] > grassBoundary)
          {
            //Generate grass
            List<Line> component = makeGrass(p, s[k]);

            lines.AddRange(component);
          }
        }
      }
    }

    //TODO GLOBAL4: here you can add extra logic to add lines that may connect between grid poitns or add other graphical elements to your design
    // Generate the path leading to visitor center
    //Start
    int startJ = minJ;
    int startI = minI;
    Point3d pStart = new Point3d(minI * dx, minJ * dy, s[minJ * rx + minI]);
    //Next
    int nextJ = minJ + 1;
    int nextI = minI;
    Point3d pNext = new Point3d(nextI * dx, nextJ * dy, s[nextJ * rx + nextI]);
    //Path and Steps
    List<Line> pathLine = new List<Line>();
    List<Line> pathStep = new List<Line>();
    //Traverse until reach the border
    while((nextJ != 0 && nextJ != (ry - 1)) && (nextI != 0 && nextI != (rx - 1)))
    {
      //traverse the 3 adjacent points(southwest, south, southeast)
      double minH = 10000.0;
      int minIndex = 0;
      for(int i = -1; i < 2; i++)
      {
        double h = Math.Abs(s[nextJ * rx + (nextI + i)] - s[startJ * rx + startI]);
        if( minH > h)
        {
          minH = h;
          minIndex = i;
        }
      }
      startJ = nextJ;
      startI = nextI + minIndex;
      pNext.X = startI * dx;
      pNext.Z = s[startJ * rx + startI];
      // Add Path Line Segment
      Line path = new Line(pStart, pNext);
      pathLine.Add(path);
      double len = 0.0;
      // Calculate the step direction
      Vector3d toEnd = (new Vector3d(pNext) - new Vector3d(pStart));
      Vector3d pNormal = Vector3d.CrossProduct(toEnd, new Vector3d(0, 0, -1.0));
      pNormal.Unitize();
      // Add Steps, step : len
      while(len < path.Length)
      {
        Point3d p = path.PointAtLength(len);
        len = len + step;
        pathStep.Add(new Line(p - pNormal * stepL, p + pNormal * stepL));
      }
      //Update Start & Next
      nextJ += 1;
      nextI = startI;
      pStart = pNext;
      pNext.Y = nextJ * dy;
    }
    //Add them to the main list
    lines.AddRange(pathLine);
    lines.AddRange(pathStep);

    // Generate Cross at the highest point on the mountain
    Point3d highest = new Point3d(maxHeightI * dx, maxHeightJ * dy, 0);
    List<Line> cross = makeCross(highest, maxHeight, focus);
    Transform scaleCross = Transform.Scale(highest + unitZ * maxHeight, 10.0);
    for (int m = 0, len = cross.Count; m < len; m++)
    {
      Line temp = cross[m];
      temp.Transform(scaleCross);
      cross[m] = temp;
    }
    lines.AddRange(cross);

    // Generate Moon:)
    int moonRadius = 10;
    double moonSegment = 0.2;
    List < Line> moon = new List<Line>();
    Random random1 = new Random(100);
    Random random2 = new Random(101);
    Random random3 = new Random(102);
    for (int i = 0; i < 1000; i++)
    {
      double moonX = random1.Next(-(moonRadius), moonRadius);
      double moonY = random2.Next(-(moonRadius), moonRadius);
      double moonZ = random3.Next(-(moonRadius), moonRadius);
      Vector3d vMoon = new Vector3d(moonX, moonY, moonZ);
      vMoon.Unitize();
      Vector3d moonNormal = Vector3d.CrossProduct(vMoon, -unitZ);
      moonNormal.Unitize();
      moon.Add(new Line(focus + vMoon + moonNormal * moonSegment, focus + vMoon - moonNormal * moonSegment));

    }
    lines.AddRange(moon);

    //END GLOBAL4

    A = lines;

  }

  // <Custom additional code> 
  //Function -> make grass//
  List<Line> makeGrass(Point3d p, double s) {
    List<Line> grass = new List<Line>();

    Point3d elevatedPoint = p;
    elevatedPoint.Z = s;

    Point3d end = elevatedPoint;
    end.Z = s + 0.2;

    grass.Add(new Line(elevatedPoint, end));

    return grass;
  }

  //Function -> make bushes//
  List<Line> makeBush(Point3d elevatedPoint, double s) {
    List<Line> bush = new List<Line>();

    Point3d top = elevatedPoint;
    top.Z += 0.5;

    Vector3d u = new Vector3d(1.0, 1.0, 0.2);
    u.Unitize();
    Point3d branch = elevatedPoint;
    branch.Z = s + 0.25;
    Point3d b1 = branch + 0.2 * u;
    bush.Add(new Line(elevatedPoint, top));
    bush.Add(new Line(branch, b1));

    Point3d start1 = top;
    Point3d start2 = b1;
    Vector3d v = new Vector3d(1.0, 1.0, 0.15);
    v.Unitize();
    for (int i = 0; i < 4; ++i) {
      Point3d end1 = start1 + 0.3 * v;
      Point3d end2 = start2 + 0.2 * v;
      bush.Add(new Line(start1, end1));
      bush.Add(new Line(start2, end2));
      start1 = end1;
      start2 = end2;
      v.Rotate(Math.PI * 0.5, top - elevatedPoint);
    }
    return bush;
  }

  //Function -> make cherry tree//
  List<Line> makeCherry(Point3d elevatedPoint, double s) {
    List<Line> cherry = new List<Line>();

    //make trunk//
    Vector3d up = new Vector3d(0.0, 0.0, 1.0);
    cherry.Add(new Line(elevatedPoint, elevatedPoint + up));

    Point3d node = new Point3d(elevatedPoint.X, elevatedPoint.Y, s + 0.3);
    Vector3d v = new Vector3d(1.0, 1.0, 1.5);
    v.Unitize();

    //make branches//
    for (int i = 0; i < 5; ++i) {
      cherry.Add(new Line(node, node + 0.5 * v));
      node.Z += 0.1;
      v.Rotate(Math.PI * 0.67, up);
      v *= 0.7;
    }

    return cherry;
  }

  //Function -> make pine tree//
  List<Line> makePine(Point3d elevatedPoint, double s) {
    List<Line> pine = new List<Line>();

    Vector3d up = new Vector3d(0.0, 0.0, 1.0);
    Vector3d node = new Vector3d(0.0, 0.0, 0.3);
    Vector3d r = new Vector3d(0.3, 0.0, 0.0);
    Vector3d f = new Vector3d(0.0, 0.3, 0.0);
    Vector3d l = -r;
    Vector3d b = -f;

    pine.Add(new Line(elevatedPoint, elevatedPoint + node));
    pine.Add(new Line(elevatedPoint + node + r, elevatedPoint + node + l));
    pine.Add(new Line(elevatedPoint + node + r, elevatedPoint + up));
    pine.Add(new Line(elevatedPoint + node + l, elevatedPoint + up));
    pine.Add(new Line(elevatedPoint + node + f, elevatedPoint + node + b));
    pine.Add(new Line(elevatedPoint + node + f, elevatedPoint + up));
    pine.Add(new Line(elevatedPoint + node + b, elevatedPoint + up));

    return pine;
  }

  //Function -> make cross//
  //let it face visitor center//
  List<Line> makeCross(Point3d p, double s, Point3d focus) {
    List<Line> cross = new List<Line>();

    Point3d elevatedPoint = p;
    elevatedPoint.Z = s;
    Point3d focusXY = focus;
    focusXY.Z = 0;
    Vector3d toFocus = focusXY - p;
    toFocus.Unitize();
    Point3d x = elevatedPoint;
    x.Z += 0.3;

    Vector3d normal = new Vector3d(0.0, 0.0, 1.0);

    Vector3d arm = toFocus;
    arm.Rotate(Math.PI * 0.5, normal);

    cross.Add(new Line(elevatedPoint, elevatedPoint + 0.4 * normal));
    cross.Add(new Line(x - 0.1 * arm, x + 0.1 * arm));

    return cross;
  }

  //Function -> decide where to put different types of plants//
  List<Line> makeComponent(Point3d p, double s, Vector3d v, Point3d focus, double maxHeight) {
    List<Line> lines = new List<Line>();
    Point3d elevatedPoint = p;
    elevatedPoint.Z = s;

    Vector3d toFocus = focus - p;

    Vector3d normal = new Vector3d(0.0, 0.0, 1.0);

    double bushBoundary = 0.35;
    double pineBoundary = 0.55;

    if (Vector3d.VectorAngle(v, normal) > 60.0) {
      //if too steep, grow grass//
      lines.AddRange(makeGrass(p, s));
    }
    else {
      if (toFocus.Length > 3 && toFocus.Length < 4.8) {
        //grow bushes around visitor center//
        lines.AddRange(makeCherry(elevatedPoint, s));
      }
      else {
        if (s > bushBoundary * maxHeight && s <= pineBoundary * maxHeight) {
          //grow bushes//
          lines.AddRange(makeBush(elevatedPoint, s));
        }
        else if (s > 0.01 * maxHeight && s <= bushBoundary * maxHeight) {
          //grow cherry trees//
          lines.AddRange(makeCherry(elevatedPoint, s));
        }
        else if (s > pineBoundary * maxHeight && s < 1.0 * maxHeight) {
          //grow pine trees//
          lines.AddRange(makePine(elevatedPoint, s));
        }
        else {
          //if too hight or too low, grow grass//
          lines.AddRange(makeGrass(p, s));
        }
      }
    }

    return lines;
  }

  // </Custom additional code> 
}
