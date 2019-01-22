using System;
using Gtk;
using OxyPlot;
using OxyPlot.GtkSharp;
using OxyPlot.Series;
using Extreme.Mathematics;
using OxyPlot.Axes;

namespace Solver
{
    class MainClass
    {
        public static void Main(string[] args)
        {
            Application.Init();
            MainWindow win = new MainWindow();
            win.Show();
            Application.Run();
        }
    }

    public partial class MainWindow : Gtk.Window
    {
        PlotView plotView;
        PlotModel plotModel;

        Extreme.Mathematics.Calculus.RombergIntegrator integrator = new Extreme.Mathematics.Calculus.RombergIntegrator();
        //SimpsonIntegrator integrator = new SimpsonIntegrator();

        double begin_domain = 0;
        double end_domain = 1;
        double k = 1;
        double u0 = 5;
        int n = 5;
        double relativeToleranceForIntegrator = 1e-6;
        int minIterationForIntegrator = 11;
        double[] linearSpace;

        public MainWindow() : base(Gtk.WindowType.Toplevel)
        {
            initialization();
            Solve();
        }

        private void initialization()
        {
            plotView = new PlotView();
            this.Add(plotView);
            plotView.ShowAll();
            plotModel = new PlotModel();

            LinearAxis linearAxisX = new LinearAxis();
            linearAxisX.MajorGridlineStyle = LineStyle.Solid;
            linearAxisX.MinorGridlineStyle = LineStyle.Dot;
            linearAxisX.Position = AxisPosition.Bottom;
            plotModel.Axes.Add(linearAxisX);

            LinearAxis linearAxisY = new LinearAxis();
            linearAxisY.MajorGridlineStyle = LineStyle.Solid;
            linearAxisY.MinorGridlineStyle = LineStyle.Dot;
            plotModel.Axes.Add(linearAxisY);

            integrator.RelativeTolerance = relativeToleranceForIntegrator;
            integrator.ConvergenceCriterion = ConvergenceCriterion.WithinRelativeTolerance;
            integrator.MinIterations = minIterationForIntegrator;

            linearSpace = new double[n + 1];
            linearSpace[0] = begin_domain;
            linearSpace[n] = end_domain;
            for(int index = 1; index < n; index++)
            {
                linearSpace[index] = (end_domain - begin_domain) / n * index;
            }
        }

        private void Solve()
        {
            Func<double, double>[] testingFunctions = new Func<double, double>[n + 1];
            var B_matrix = MathNet.Numerics.LinearAlgebra.Matrix<double>.Build.Dense(n, n);
            var B_matrix2 = MathNet.Numerics.LinearAlgebra.Matrix<double>.Build.Dense(n, n);
            var L_matrix = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(n);

            for (int seria = 0; seria <= n; seria++)
            {
                int uniqSeriaCopy = seria;
                testingFunctions[seria] = xVar => Math.Max(0, 1 - (Math.Abs(xVar - (end_domain - begin_domain) * uniqSeriaCopy / n) * n));

                plotModel.Series.Add(new FunctionSeries(testingFunctions[seria], 0, 1, 500));
            }

            for (int firstIndex = 0; firstIndex < n; firstIndex++)
            {
                for (int secondIndex = 0; secondIndex < n; secondIndex++)
                {
                    B_matrix[firstIndex, secondIndex] = B_function(testingFunctions[secondIndex + 1], testingFunctions[firstIndex + 1]);
                }
            }

            for (int index = 0; index < n; index++)
            {
                L_matrix[index] = L_function(testingFunctions[index + 1]) - u0 * B_function(testingFunctions[0], testingFunctions[index + 1]);
            }

            var U_matrix = B_matrix.Solve(L_matrix);
            var new_U_matrix = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(n + 1);
            new_U_matrix[0] = u0;
            for (int i = 1; i < n + 1; i++)
            {
                new_U_matrix[i] = U_matrix[i - 1];
            }
            /*
            var U_matrix2 = B_matrix2.Solve(L_matrix);
            var new_U_matrix2 = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(n + 1);
            new_U_matrix2[0] = u0;
            for (int i = 1; i < n + 1; i++)
            {
                new_U_matrix2[i] = U_matrix2[i - 1];
            }*/

            Console.WriteLine("B matrix");
            Console.WriteLine(B_matrix);
            Console.WriteLine("L matrix");
            Console.WriteLine(L_matrix);
            Console.WriteLine("U matrix");
            Console.WriteLine(U_matrix);

            /*Console.WriteLine("B matrix");
            Console.WriteLine(B_matrix2);
            Console.WriteLine("L matrix");
            Console.WriteLine(L_matrix);
            Console.WriteLine("U matrix");
            Console.WriteLine(U_matrix);*/

            LineSeries approximFunction = new LineSeries();
            //LineSeries approximFunction2 = new LineSeries();

            foreach (double point in linearSpace)
            {
                double valInX = 0;
                //double valInX2 = 0;
                for (int index = 0; index < n + 1; index++)
                {
                    valInX += new_U_matrix[index] * testingFunctions[index](point);
                    //valInX2 += new_U_matrix2[index] * testingFunctions[index](point);
                }
                approximFunction.Points.Add(new DataPoint(point, valInX));
                //approximFunction2.Points.Add(new DataPoint(point, valInX2));
            }

            plotModel.Series.Add(approximFunction);
            //plotModel.Series.Add(approximFunction2);
            plotView.Model = plotModel;
        }

        private double B_function(Func<double, double> u, Func<double, double> v)
        {
            Func<double, double> dU = u.GetNumericalDifferentiator();
            Func<double, double> dV = v.GetNumericalDifferentiator();

            Func<double, double> partial1 = x => dU(x) * dV(x);
            Func<double, double> partial2 = x => dU(x) * v(x);

            //k * dU(0) * v(0)
            var result = k * integrator.Integrate(partial1, begin_domain, end_domain)
                + integrator.Integrate(partial2, begin_domain, end_domain);
                
            return result;
        }

        private double L_function(Func<double, double> v)
        {
            Func<double, double> function = x => (5 * x - 10) * v(x); 

            var result = integrator.Integrate(function, begin_domain, end_domain) + 3 * k * v(1);
            return result;
        }

        protected void OnDeleteEvent(object sender, DeleteEventArgs a)
        {
            Application.Quit();
            a.RetVal = true;
        }

    }
}
