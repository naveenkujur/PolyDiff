using Nori;

/// <summary>An abstract alignment pair of spans along PolyA and PolyB</summary>
/// <param name="StartLieA">Start lie along Poly A</param>
/// <param name="EndLieA">End lie along Poly A</param>
/// <param name="StartLieB">Start lie along Poly B</param>
/// <param name="EndLieB">End lie along Poly B</param>
/// <remarks>
/// Used to model both types of alignment span pairs:
///   the persistent (overlapping) span pair
///   and the substitution (replacement) span pair
/// </remarks>
readonly struct PairedSpan { // Persistent data | But could be initialized using SegSlice/s
   public PairedSpan (Seg segA, double startLieA, double endLieA, Seg segB, double startLieB, double endLieB)
      => (SegA, StartLieA, EndLieA, SegB, StartLieB, EndLieB)
      = (segA, (float)startLieA, (float)endLieA, segB, (float)startLieB, (float)endLieB);

   readonly Seg SegA, SegB;
   public readonly float StartLieA, EndLieA;
   public readonly float StartLieB, EndLieB;

   public SegSlice SegSliceA => new (SegA, StartLieA, EndLieA);
   public SegSlice SegSliceB => new (SegB, StartLieB, EndLieB);

   public override string ToString () => $"[({StartLieA},{EndLieA})({StartLieB},{EndLieB})]";
}

class Coverage {
   Coverage (Poly a, Poly b) { }
   // Coverage == List of overlapping spans
   //    This is sufficient for overlap cleanup operation
   // Difference b/w two polys == PolyB - PolyA
   //    Shows how PolyA changes to PolyB
   //       How to go from PolyA to Polyb
   //          Which spans persist and get carried over, and which spans undergo span substitution!
   //    There may be overlapping spans which remain - coverage - persistent spans
   //    There are spans in PolyA which fade in - outgoing spans
   //    There are spans in PolyB which fade out - incoming spans
   public static Coverage? Compute (Poly polyA, Poly polyB) {
      if ((polyA.GetBound ().InflatedL (Lib.Epsilon) * polyB.GetBound ()).IsEmpty) return null;

      // NOTE: All Polys are essentially treated like open polys in this case of overlap detection and ranging.
      var cvrg = new Coverage (polyA, polyB);
      SegSlice sa = new (polyA[0]), sb = new (polyB[0]);
      PairedSpan? overlap0 = null; // Hold the first overlap; This will help determine the STOP condition.
      while (true) {
         if (overlap0 is { } o0) {
            // Slice already part of first overlap, we're done!
            if (o0.SegSliceA.EQ (sa)) break;
         }
         var overlap = sa.ComputeOverlap (sb);
         if (overlap is { } o) {
            if (overlap0 is null) overlap0 = o; // Register the first overlap
            cvrg.mOverlaps.Add (o);
            // Get the slices from overlap!
            sa = o.SegSliceA;
            sb = o.SegSliceB;
            continue;
         }
         if (sb.IsLast) {
            // PolyB is consumed.
            // Do we reenter PolyB
            if (overlap0 is null) { // No overlap has yet been detected
               if (sa.IsLast) break; // No overlap b/w PolyA and PolyB
               // Continue with the next seg and look for overlap with PolyB again.
               sa = sa.Next ();
               sb = sb.First (); // Reset
               continue;
            }
            // Some overlap has been found
            sb = sb.First (); // Reset
         } else sb = sb.Next ();
         if (sa.IsLast) {
            if (overlap0 is null) break;
            sa = sa.First (); // Reset
         }
      }
      return cvrg;
   }

   public IReadOnlyList<PairedSpan> Overlaps => mOverlaps;

   readonly List<PairedSpan> mOverlaps = [];
}

readonly ref struct SegSlice {
   public SegSlice () => throw new NotImplementedException ();
   /// <summary>Ctor</summary>
   /// <param name="seg">Sliced seg</param>
   /// <param name="startLie">Slice start lie, along the underlying seg</param>
   /// <param name="endLie">Slice end lie, along the underlying seg</param>
   public SegSlice (Seg seg, double startLie = 0, double endLie = 1) => (mSeg, mStartLie, mEndLie) = (seg, startLie, endLie);

   readonly Seg mSeg;
   readonly double mStartLie, mEndLie;

   public bool IsLast => mSeg.IsLast && mEndLie.EQ (1);
   public SegSlice Next () => mEndLie < (1 - Lib.Epsilon) ? new SegSlice (mSeg, mEndLie, 1)
      : new SegSlice (mSeg.Next!.Value);
   public SegSlice First () => new (mSeg.Poly[0]);

   public bool EQ (SegSlice o)
      => mSeg.Poly == o.mSeg.Poly && mSeg.N == o.mSeg.N && mStartLie.EQ (o.mStartLie); // No need to match endLie!

   public PairedSpan? ComputeOverlap (SegSlice o) {
      if (!IsAlignedWith (o)) return null;
      var (sa, sb) = (mSeg, o.mSeg);
      if (mSeg.IsLine) {
         var (sLie, eLie) = (sa.GetLie (sb.A), sa.GetLie (sb.B)); // Lie of segB, in terms of segA
         if (sLie > (1 - Lib.Epsilon) || eLie < Lib.Epsilon) return null;
         // Got overlap!
         var (sLie2, eLie2) = (sb.GetLie (sa.A), sb.GetLie (sa.B)); // Lie of segA, in terms of segB
         return new PairedSpan (sa, sLie, eLie, sb, sLie2, eLie2);
      }
      {
         // Note: Arc external lies are tricky!
         // That is, if the lie is within [0..1] we're good. Otherwise, we are in a very deep soup!

         var (sLie, eLie) = (sa.GetLie (sb.A), sa.GetLie (sb.B)); // Lie of segB, in terms of segA
         if (sLie > (1 - Lib.Epsilon) || eLie < Lib.Epsilon) return null;

         // Clamp the projected lies to 0..1 range, very carefully!
         if (!InRange (sLie)) sLie = 0;
         if (!InRange (eLie)) eLie = 1;


         var (sLie2, eLie2) = (sb.GetLie (sa.A), sb.GetLie (sa.B)); // Lie of segA, in terms of segB
         // Clamp the projected lies to 0..1 range, very carefully!
         if (!InRange (sLie2)) sLie2 = 0;
         if (!InRange (eLie2)) eLie2 = 1;

         return new PairedSpan (sa, sLie, eLie, sb, sLie2, eLie2);

         static bool InRange (double lie) => lie > Lib.Epsilon && lie < 1 - Lib.Epsilon;
      }
   }


   // Implementation
   bool IsAlignedWith (SegSlice o) { // Alignment check: Common line def | common circle def
      if (mSeg.IsLine ^ o.mSeg.IsArc) return false;
      var (sa, sb) = (mSeg, o.mSeg);
      if (mSeg.IsLine) {
         if (sa.A.Side (sb.A, sb.B) != 0 || sa.B.Side (sb.A, sb.B) != 0) return false;
         if (!sa.Slope.EQ (sb.Slope)) throw new NotSupportedException ("Same direction only, for now");
         return true;
      }
      if (!sa.Center.EQ (sb.Center) || !sa.Radius.EQ (sb.Radius)) return false;
      if (sa.IsCCW != sb.IsCCW) throw new NotSupportedException ("Same winding only, for now");
      return true;
   }
}

static class Q {
   public static void RunTests () {
      TestLineSeg ();
   }

   // Single line segment modified
   public static void TestLineSeg () {
      Poly a = Poly.Parse ("M0,0 H100");
      {
         Poly b = a; // Full overlap!
         Validate (a, b, [new (a[0], 0, 1, b[0], 0, 1)]);
      }
      {
         Poly b = Poly.Parse ("M0,0 H200"); // Elongated line
         Validate (a, b, [new (a[0], 0, 1, b[0], 0, 0.5)]);
      }
      {
         Poly b = Poly.Parse ("M10,0 H110"); // Shifted forward
         Validate (a, b, [new (a[0], 0.1, 1, b[0], 0, 0.9)]);
      }
      {
         Poly b = Poly.Parse ("M-10,0H90"); // Shifted backward
         Validate (a, b, [new (a[0], 0, 0.9, b[0], 0.1, 1)]);
      }
      {
         Poly b = Poly.Parse ("M10,0H90"); // Fully contained
         Validate (a, b, [new (a[0], 0.1, 0.9, b[0], 0, 1)]);
      }
      {
         Poly a2 = Poly.Parse ("M10,0H90");
         Poly b = Poly.Parse ("M0,0H100"); // Fully contains other
         Validate (a2, b, [new (a[0], 0, 1, b[0], 0.1, 0.9)]);
      }
   }

   // Performs regression testing
   public static void Validate (Poly a, Poly b, IList<PairedSpan>? refRanges = null) {
      var cvrg = Coverage.Compute (a, b);
      var ranges = cvrg!.Overlaps;
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
      foreach (var range in ranges)
         Console.Write (range);
      Console.WriteLine ($"Pass");
   }
}
