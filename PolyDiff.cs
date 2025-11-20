using Nori;

/// <summary>Overlapping ranges along Poly A and B</summary>
/// <param name="StartLieA">Range start along Poly A</param>
/// <param name="EndLieA">Range end along Poly A</param>
/// <param name="StartLieB">Range start along Poly B</param>
/// <param name="EndLieB">Range start along Poly B</param>
readonly struct ORange {
   public ORange (double startLieA, double endLieA, double startLieB, double endLieB)
      => (StartLieA, EndLieA, StartLieB, EndLieB)
      = (startLieA.Clamp (), endLieA.Clamp (), startLieB.Clamp (), endLieB.Clamp ());

   public readonly double StartLieA;
   public readonly double EndLieA;
   public readonly double StartLieB;
   public readonly double EndLieB;

   public override string ToString () => $"[({StartLieA},{EndLieA})({StartLieB},{EndLieB})]";
   public ORange Rebased (int segAIdx, int segBIdx)
      => new (segAIdx + StartLieA, segAIdx + EndLieA, segBIdx + StartLieB, segBIdx + EndLieB);
}

static class Q {
   public static IList<ORange> PolyDiff (Poly a, Poly b) {
      if (!a.IsOpen || !b.IsOpen) throw new NotSupportedException ("Open polys only, for now");
      List<ORange> overlaps = [];
      for (int i = 0; i < a.Count; i++) {
         var segA = a[i];
         for (int j = 0; j < b.Count; j++) {
            var segB = b[j];
            if (segA.IsLine != segB.IsLine) { continue; }
            // Either both are lines, or both are arcs
            if (segA.IsLine) {
               var range = CollinearOverlap (segA, segB);
               if (range == null) continue;
               overlaps.Add (range.Value.Rebased (segA.N, segB.N));
            } 
         }
      }
      return overlaps;
   }

   /// <summary>Collinear overlap range computation</summary>
   public static ORange? CollinearOverlap (Seg segA, Seg segB) {
      if (!segA.IsLine || !segB.IsLine) throw new InvalidOperationException ();
      if (segA.A.Side (segB.A, segB.B) != 0 || segA.B.Side (segB.A, segB.B) != 0) return null; 
      if (!segA.Slope.EQ (segB.Slope)) throw new NotSupportedException ("Same direction only, for now");
      // Check for overlap
      var (sLie, eLie) = (segA.GetLie (segB.A), segA.GetLie (segB.B)); // Lie of segB, in terms of segA
      if (sLie > 1 - Lib.Epsilon || eLie < Lib.Epsilon) return null;
      // Got overlap!
      var (sLie2, eLie2) = (segB.GetLie (segA.A), segB.GetLie (segA.B)); // Lie of segA, in terms of segB
      return new ORange (sLie, eLie, sLie2, eLie2);
   }

   /// <summary>Concentric overlap range computation</summary>
   public static ORange? ConcentricOverlap (Seg segA, Seg segB) {
      if (segA.IsCircle || segB.IsCircle) throw new NotSupportedException ("No circles, for now");
      if (!segA.IsArc || !segB.IsArc) throw new InvalidOperationException ();
      if (!segA.Center.EQ (segB.Center) || !segA.Radius.EQ (segB.Radius)) return null;
      if (segA.IsCCW != segB.IsCCW) throw new NotSupportedException ("Same winding only, for now");

      var (sLie, eLie) = (segA.GetLie (segB.A), segA.GetLie (segB.B)); // Lie of segB, in terms of segA
      
      // Note: Arc external lies are tricky!
      // So we simply check if any of the projected lies are in 0..1 range!
      bool overlap = InRange (sLie) || InRange (eLie);
      if (!overlap) return null;

      // Clamp the projected lies to 0..1 range, very carefully!
      if (!InRange (sLie)) sLie = 0;
      if (!InRange (eLie)) eLie = 1;

      var (sLie2, eLie2) = (segB.GetLie (segA.A), segB.GetLie (segA.B)); // Lie of segA, in terms of segB
      // Clamp the projected lies to 0..1 range, very carefully!
      if (!InRange (sLie2)) sLie2 = 0;
      if (!InRange (eLie2)) eLie2 = 1;

      return new ORange (sLie, eLie, sLie2, eLie2);

      static bool InRange (double lie) => lie > Lib.Epsilon && lie < 1 - Lib.Epsilon;
   }

   public static void RunTests () {
      TestLineSeg ();
   }

   // Single line segment modified
   public static void TestLineSeg () {
      Poly a = Poly.Parse ("M0,0 H100");
      {
         Poly b = a; // Full overlap!
         Validate (a, b, [new (0, 1, 0, 1)]);
      }
      {
         Poly b = Poly.Parse ("M0,0 H200"); // Elongated line
         Validate (a, b, [new (0, 1, 0, 0.5)]);
      }
      {
         Poly b = Poly.Parse ("M10,0 H110"); // Shifted forward
         Validate (a, b, [new (0.1, 1, 0, 0.9)]);
      }
      {
         Poly b = Poly.Parse ("M-10,0H90"); // Shifted backward
         Validate (a, b, [new (0,0.9,0.1,1)]);
      }
      {
         Poly b = Poly.Parse ("M10,0H90"); // Fully contained
         Validate (a, b, [new (0.1,0.9,0,1)]);
      }
      {
         Poly a2 = Poly.Parse ("M10,0H90");
         Poly b = Poly.Parse ("M0,0H100"); // Fully contains other
         Validate (a2, b, [new (0,1,0.1,0.9)]);
      }
   }

   // Performs regression testing
   public static void Validate (Poly a, Poly b, IList<ORange>? refRanges = null) {
      var ranges = PolyDiff (a, b);
      if (refRanges == null) {
         foreach (var range in ranges)
            Console.WriteLine (range.ToString ());
         return;
      }
      foreach (var t in Enumerable.Zip (refRanges, ranges)) {
         if (!t.First.StartLieA.EQ (t.Second.StartLieA) || !t.First.EndLieA.EQ (t.Second.EndLieA)
            || !t.First.StartLieB.EQ (t.Second.StartLieB) || !t.First.EndLieB.EQ (t.Second.EndLieB)) {
            Console.WriteLine ($"Fail: {t}");
            return;
         }
      }
      foreach(var range in ranges) 
         Console.Write (range);
      Console.WriteLine ($"Pass");
   }
}
