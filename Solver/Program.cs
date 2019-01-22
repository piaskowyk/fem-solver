using System;
using System.IO;
using Extreme.Mathematics;
using Extreme.Mathematics.Calculus;
using Gtk;
using Newtonsoft.Json;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.GtkSharp;
using OxyPlot.Series;

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

        double begin_domain = 0;
        double end_domain = 1;
        double u0 = 5;
        double[] linearSpace;
        Config config;

        public MainWindow() : base(Gtk.WindowType.Toplevel)
        {
            initialization();
            Solve();
        }

        private void initialization()
        {
            string jsonConf = File.ReadAllText("config.json");
            config = JsonConvert.DeserializeObject<Config>(jsonConf);

            plotView = new PlotView();
            this.Add(plotView);
            plotView.ShowAll();
            plotModel = new PlotModel();

            LinearAxis linearAxisX = new LinearAxis
            {
                MajorGridlineStyle = LineStyle.Solid,
                MinorGridlineStyle = LineStyle.Dot,
                Position = AxisPosition.Bottom
            };
            plotModel.Axes.Add(linearAxisX);

            LinearAxis linearAxisY = new LinearAxis
            {
                MajorGridlineStyle = LineStyle.Solid,
                MinorGridlineStyle = LineStyle.Dot
            };
            plotModel.Axes.Add(linearAxisY);

            integrator.RelativeTolerance = config.relativeToleranceForIntegrator;
            integrator.ConvergenceCriterion = ConvergenceCriterion.WithinRelativeTolerance;
            integrator.MinIterations = config.minIterationForIntegrator;

            linearSpace = new double[config.n + 1];
            linearSpace[0] = begin_domain;
            linearSpace[config.n] = end_domain;
            for(int index = 1; index < config.n; index++)
            {
                linearSpace[index] = (end_domain - begin_domain) / config.n * index;
            }
        }

        private void Solve()
        {
            Func<double, double>[] basicFunctions = new Func<double, double>[config.n + 1];
            var B_matrix = MathNet.Numerics.LinearAlgebra.Matrix<double>.Build.Dense(config.n, config.n);
            var L_matrix = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(config.n);

            for (int seria = 0; seria <= config.n; seria++)
            {
                int uniqSeriaCopy = seria;
                basicFunctions[seria] = xVar => 
                    Math.Max(0, 1 - (Math.Abs(xVar - (end_domain - begin_domain) * uniqSeriaCopy / config.n) * config.n));

                if (config.showBasicFunctions)
                {
                    plotModel.Series.Add(new FunctionSeries(basicFunctions[seria], begin_domain, end_domain, config.n * 100));
                }
            }

            for (int firstIndex = 0; firstIndex < config.n; firstIndex++)
            {
                for (int secondIndex = 0; secondIndex < config.n; secondIndex++)
                {
                    B_matrix[firstIndex, secondIndex] = B_function(basicFunctions[secondIndex + 1], basicFunctions[firstIndex + 1]);
                }
            }

            for (int index = 0; index < config.n; index++)
            {
                L_matrix[index] = L_function(basicFunctions[index + 1]) - u0 * B_function(basicFunctions[0], basicFunctions[index + 1]);
            }

            var U_matrix = B_matrix.Solve(L_matrix);
            var increased_U_matrix = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(config.n + 1);
            increased_U_matrix[0] = u0;
            for (int i = 1; i < config.n + 1; i++)
            {
                increased_U_matrix[i] = U_matrix[i - 1];
            }

            Console.WriteLine("B matrix");
            Console.WriteLine(B_matrix);
            Console.WriteLine("L matrix");
            Console.WriteLine(L_matrix);
            Console.WriteLine("U matrix");
            Console.WriteLine(U_matrix);

            LineSeries approximFunction = new LineSeries();

            foreach (double point in linearSpace)
            {
                double valInX = 0;
                for (int index = 0; index < config.n + 1; index++)
                {
                    valInX += increased_U_matrix[index] * basicFunctions[index](point);
                }
                approximFunction.Points.Add(new DataPoint(point, valInX));
            }

            plotModel.Series.Add(approximFunction);
            plotView.Model = plotModel;
        }

        private double B_function(Func<double, double> u, Func<double, double> v)
        {
            Func<double, double> dU = u.GetNumericalDifferentiator();
            Func<double, double> dV = v.GetNumericalDifferentiator();

            double partial1(double x) => dU(x) * dV(x);
            double partial2(double x) => dU(x) * v(x);

            var result = config.k * integrator.Integrate(partial1, begin_domain, end_domain)
                + integrator.Integrate(partial2, begin_domain, end_domain);
                
            return result;
        }

        private double L_function(Func<double, double> v)
        {
            double function(double x) => (5 * x - 10) * v(x);

            var result = integrator.Integrate(function, begin_domain, end_domain) + 3 * config.k * v(1);
            return result;
        }

        protected void OnDeleteEvent(object sender, DeleteEventArgs a)
        {
            Application.Quit();
            a.RetVal = true;
        }

        public class Config
        {
            public int k;
            public int n;
            public double relativeToleranceForIntegrator;
            public int minIterationForIntegrator;
            public bool showBasicFunctions;
        }

    }
}
