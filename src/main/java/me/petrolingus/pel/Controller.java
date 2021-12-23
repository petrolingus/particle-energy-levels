package me.petrolingus.pel;

import javafx.scene.chart.LineChart;
import javafx.scene.chart.XYChart;
import javafx.scene.control.TextField;

import java.util.ArrayList;
import java.util.List;

public class Controller {

    public TextField sigmaText;
    public TextField pitWidthText;
    public TextField kText;
    public TextField energyText;
    public TextField kthText;

    double U0;
    double koef;
    int k_dlg;
    double R;
    double a;
    double Emin;
    double Emax;

    double[] psi;
    double[] Fi;
    double[] r;
    double[] lastFi;
    double Sigma;

    public LineChart<Number, Number> phaseEnergyChart;
    public LineChart<Number, Number> waveFunctionChart;

    double calculateU(double koef, double x) {
        double c0 = U0 * Math.exp(-(x * x) / (2 * Sigma * Sigma));
        double c1 = Sigma * Math.sqrt(2 * Math.PI);
//        return 0.5 * koef * x * x + c0 / c1;
        return c0 / c1;
    }

    double dFi(double energy, double fi, double z) {
        return (calculateU(koef, z) - energy) * Math.pow(Math.cos(fi), 2) - Math.pow(Math.sin(fi), 2);
    }

    double getFi(double energy, int count) {
        double[] mas = new double[count];
        double k1, k2, k3, k4;
        double zStart = -R;
        double zFinal = R;

        int CurIndex = 0;
        double CurZ = 0;
        double dz = (zFinal - zStart) / count;
        mas[CurIndex] = Math.PI / 2;
        Fi[CurIndex] = mas[CurIndex];

        for (int i = 1; i < count; i++) {
            CurZ = zStart + i * dz;

            k1 = dFi(energy, mas[CurIndex], CurZ) * dz;
            k2 = dFi(energy, mas[CurIndex] + k1 / 2, CurZ + dz / 2) * dz;
            k3 = dFi(energy, mas[CurIndex] + k2 / 2, CurZ + dz / 2) * dz;
            k4 = dFi(energy, mas[CurIndex] + k3, CurZ + dz) * dz;

            mas[CurIndex + 1] = mas[CurIndex] + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
            Fi[CurIndex] = mas[CurIndex];

            CurIndex++;
        }

        return mas[count - 1];
    }

    double getEstat(double eMin, double eMax, double fiLeft, double fiRight, double statValue, double eps, int count) {
        double midE;
        double midFi;
        int iter = 0;
        double err;

        do {
            midE = (eMax + eMin) / 2;
            midFi = getFi(midE, count);
            double flag = (fiLeft - statValue) * (midFi - statValue);

            if (flag < 0) {
                eMax = midE;
                fiRight = midFi;
            } else {
                eMin = midE;
                fiLeft = midFi;
            }
            err = Math.abs(midFi - statValue);
        }
        while (err > eps);

        return midE;
    }

    double[] getR(double _E, int count) {
        double[] mas = new double[count];

        double zStart = -R;
        double zFinal = R;

        int CurIndex = count / 2;
        double CurZ;
        double dz = (zFinal - zStart) / count;
        double k1, k2, k3, k4;
        mas[CurIndex] = 1;

        for (int i = CurIndex + 1; i < count - 1; i++) {
            CurZ = zStart + i * dz;
            getFi(_E, count);
            k1 = dz * dr(mas[CurIndex], calculateU(koef, CurZ), _E, Fi[CurIndex]);
            k2 = dz * dr(mas[CurIndex] + k1 / 2, calculateU(koef, CurZ + dz / 2), _E, Fi[CurIndex]);
            k3 = dz * dr(mas[CurIndex] + k2 / 2, calculateU(koef, CurZ + dz / 2), _E, Fi[CurIndex]);
            k4 = dz * dr(mas[CurIndex] + k3, calculateU(koef, CurZ + dz), _E, Fi[CurIndex]);

            mas[CurIndex + 1] = mas[CurIndex] + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
            CurIndex++;
        }

        dz = -(zFinal - zStart) / count;
        CurIndex = count / 2 - 1;
        mas[CurIndex] = 1;

        for (int i = CurIndex - 1; i >= 0; i--) {
            CurZ = zStart - i * dz;
            getFi(_E, count);
            k1 = dr(mas[CurIndex], calculateU(koef, CurZ), _E, Fi[CurIndex]);
            k2 = dr(mas[CurIndex] + k1 * dz / 2, calculateU(koef, CurZ + dz / 2), _E, Fi[CurIndex]);
            k3 = dr(mas[CurIndex] + k2 * dz / 2, calculateU(koef, CurZ + dz / 2), _E, Fi[CurIndex]);
            k4 = dr(mas[CurIndex] + k3 * dz, calculateU(koef, CurZ + dz), _E, Fi[CurIndex]);

            CurIndex--;
            mas[CurIndex] = mas[CurIndex + 1] + (k1 + 2 * k2 + 2 * k3 + k4) * dz / 6;

        }
        return mas;
    }

    double dr(double r, double U, double E, double fi) {
        return r * (U - E + 1) * Math.cos(fi) * Math.sin(fi);
    }

    public void onDrawButton() {

        koef = Double.parseDouble(kText.getText());
        a = Double.parseDouble(pitWidthText.getText());
        R = 15;
        U0 = Double.parseDouble(energyText.getText());
        Sigma = Double.parseDouble(sigmaText.getText());

        int count = 300;
        psi = new double[count];
        Fi = new double[count];
        r = new double[count];
        lastFi = new double[count];

        Emin = 0;
        Emax = koef * a;

        double currentFi = -Math.PI / 2;
        double currentE = Emin;

        lastFi[0] = getFi(currentE, count);

        List<Double> energyLevels = new ArrayList<>();

        int k = 0;

        for (int i = 1; i < count; i++) {
            currentE = Emin + ((Emax - Emin) / count) * i;
            lastFi[i] = getFi(currentE, count);
            if (lastFi[i] < currentFi && lastFi[i - 1] > currentFi) {
                k++;
                double energy = getEstat(currentE - (Emax - Emin) / count, currentE, lastFi[i - 1], lastFi[i], currentFi, 0.1, count);
                energyLevels.add(energy);
                currentFi = -(2 * k + 1) * Math.PI / 2;
            }
        }

        System.out.println(k);

        // Draw phase chart
        XYChart.Series<Number, Number> series1 = new XYChart.Series<>();
        for (int i = 0; i < count; i++) {
            series1.getData().add(new XYChart.Data<>(i, lastFi[i]));
        }
        phaseEnergyChart.getData().clear();
        phaseEnergyChart.getData().add(series1);

        // Draw wave function chart
        waveFunctionChart.getData().clear();
        for (Double energyLevel : energyLevels) {
            XYChart.Series<Number, Number> waveFunction = new XYChart.Series<>();
            r = getR(energyLevel, count);
            getFi(energyLevel, count);
            for (int j = 0; j < count; j++) {
                double value = r[j] * Math.cos(Fi[j]);
                waveFunction.getData().add(new XYChart.Data<>(j, value));
            }
            waveFunctionChart.getData().add(waveFunction);
        }
    }
}
