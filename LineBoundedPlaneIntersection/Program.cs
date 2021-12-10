using MathNet.Numerics.LinearAlgebra;
using System;
using System.Numerics;

namespace LineBoundedPlaneIntersection
{
    class Program
    {
        static void Main(string[] args)
        {
            Vector3 u = new(1, 0, 0);
            Vector3 v = new(0, 1, 0);
            Plane plane = new(new(0, 0, 0.25f), u, v);
         
            Vector3 a = new(0.5f, 0.5f, 0.5f);
            Vector3 b = new(0.5f, 0.5f, 1.0f);
            LineSegment line = new(a, b);

            IntersectHit hit = Intersect.FindPoint(plane, line);

            Console.WriteLine(hit);


            Triangle triangle = new(new(0, 0, 0), new(1, 0, 0), new(0, 1, 0));
            IntersectHit tHit = Intersect.FindPoint(triangle, line);

            Console.WriteLine(tHit);
        }
    }

    static class Intersect
    {
        // Line/Plane intersection: https://stackoverflow.com/a/18543221/2223566
        public static IntersectHit FindPoint(Plane plane, LineSegment line)
        {
            Vector3 U = line.B - line.A;
            Vector3 PN = plane.GetNormal();
            float dot = Vector3.Dot(PN, U);
            float eps = 1e-6f;

            if (Math.Abs(dot) > eps)
            {
                Vector3 W = line.A - plane.Position;
                float fac = -Vector3.Dot(PN, W) / dot;
                Vector3 P = line.A + U * fac;

                if(!line.IsPointOnLineSegment(P))
                {
                    return new(Vector3.Zero, IsHit: false, Reason: "Point lies outside of the line segment.");
                }

                // Check if the point can be expressed as a linear combination of the span of the plane, whilst
                // the scalars are both in the range [0, 1].
                if(!plane.IsPointInFinitePlane(P))
                {
                    return new(Vector3.Zero, IsHit: false, Reason: "Point lies outside of finite plane.");
                }


                return new(P, IsHit: true);
                
            }

            return new(Vector3.Zero, IsHit: false, Reason: "Line is parallel to the plane.");
        }

        // https://stackoverflow.com/a/42752998
        public static IntersectHit FindPoint(Triangle triangle, LineSegment line)
        {
            Vector3 RD = Vector3.Normalize(line.B - line.A);

            Vector3 A = triangle.A;
            Vector3 B = triangle.B;
            Vector3 C = triangle.C;

            Vector3 E1 = B - A;
            Vector3 E2 = C - A;
            Vector3 N = Vector3.Cross(E1, E2);
            float det = -Vector3.Dot(RD, N);
            float invdet = 1.0f / det;
            Vector3 AO = line.A - A;
            Vector3 DAO = Vector3.Cross(AO, RD);
            float u = Vector3.Dot(E2, DAO) * invdet;
            float v = -Vector3.Dot(E1, DAO) * invdet;
            float t = Vector3.Dot(AO, N) * invdet;
            float eps = 1e-6f;
            bool isHit = Math.Abs(det) >= eps && u >= 0.0 && v >= 0.0 && (u + v) <= 1.0;

            if(!isHit)
            {
                return new(Vector3.Zero, IsHit: false, Reason: "No hit.");
            }

            Vector3 P = line.A + t * RD;

            if (!line.IsPointOnLineSegment(P))
            {
                return new(Vector3.Zero, IsHit: false, Reason: "Point lies outside of the line segment.");
            }

            return new(P, IsHit: true);
        }
    }

    record IntersectHit(Vector3 Point, bool IsHit, string Reason = "");

    class Triangle
    {
        public Vector3 A { get; set; }
        public Vector3 B { get; set; }
        public Vector3 C { get; set; }
        public Triangle(Vector3 a, Vector3 b, Vector3 c)
        {
            A = a;
            B = b;
            C = c;
        }
    }

    class Plane
    {
        public Vector3 Position { get; }
        public Vector3 U { get; }
        public Vector3 V { get; }

        public Plane(Vector3 position, Vector3 u, Vector3 v)
        {
            Position = position;
            U = u;
            V = v;
        }

        public Vector3 GetPoint(float u, float v)
        {
            return Position + u * U + v * V;
        }

        public Vector3 GetNormal()
        {
            return Vector3.Cross(U, V);
        }

        // Check if we can make a linear combination of U and V, where alpha en beta are the scalars.
        // If the scalars are both between [0, 1] we are inside the finite plane.
        public bool IsPointInFinitePlane(Vector3 P)
        {
            // Solve the equation:  P = alpha * v1 + beta * v2, where alpha and beta are unknown.
            var M = Matrix<float>.Build.DenseOfArray(new float[,] {
                { U.X, V.X },
                { U.Y, V.Y },
                { U.Z, V.Z }
            });
            
            var b = MathNet.Numerics.LinearAlgebra.Vector<float>.Build.Dense(new float[] { P.X, P.Y, P.Z });
            var x = M.Solve(b);

            float alpha = x[0];
            float beta = x[1];

            float eps = 1e-6f;

            if(alpha - eps < 0 || alpha + eps > 1)
            {
                return false;
            }

            if(beta - eps < 0 || beta + eps > 1)
            {
                return false;
            }

            return true;
        }
    }

    class LineSegment
    {
        public Vector3 A { get; }
        public Vector3 B { get; }

        public LineSegment(Vector3 a, Vector3 b)
        {
            A = a;
            B = b;
        }

        public Vector3 GetPoint(float scalar)
        {
            var direction = B - A;
            return A + scalar * direction;
        }

        // https://lucidar.me/en/mathematics/check-if-a-point-belongs-on-a-line-segment/
        public bool IsPointOnLineSegment(Vector3 C)
        {
            float eps = 1e-6f;

            // Check if the vectors are colinear to each other, the cross product will
            // return zero in this case.
            Vector3 cross = Vector3.Cross(B - A, C - A);

            // If AC and BC are colinear, then the cross product ABxAC=0, if not
            // then the point is not on the line.
            if(Vector3.Dot(cross, cross) >= eps)
            {
                return false;
            }

            float K_AC = Vector3.Dot(B - A, C - A);
            float K_AB = Vector3.Dot(B - A, B - A);

            // Point is not between A and B.
            if(K_AC + eps < 0)
            {
                return false;
            }

            // Point is not between A and B.
            if(K_AC > K_AB + eps)
            {
                return false;
            }

            // Point is on A, or B, or inbetween AB.
            return K_AC + eps > 0 && K_AC < K_AB + eps;
        }
    }
}
