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

    private static final int COUNT = 300; // Count of samples
    private static final double R = 15; // Z-interval
    private static final double EPS = 0.1; // Z-interval

    double U0;
    double coefficient;
    double a;
    double sigma;

    double[] fi;
    double[] r;
    double[] lastFi;

    public LineChart<Number, Number> phaseEnergyChart;
    public LineChart<Number, Number> waveFunctionChart;

    double calculateU(double a, double x) {
        double c0 = U0 * Math.exp(-(x * x) / (2 * sigma * sigma));
        double c1 = sigma * Math.sqrt(2 * Math.PI);
        return 0.5 * a * x * x + c0 / c1;
    }

    double getFi(double energy) {
        double[] mas = new double[COUNT];
        double k1, k2, k3, k4;
        double zStart = -R;

        int CurIndex = 0;
        double CurZ;
        double dz = (R - zStart) / COUNT;
        mas[CurIndex] = Math.PI / 2;
        fi[CurIndex] = mas[CurIndex];

        for (int i = 1; i < COUNT; i++) {
            CurZ = zStart + i * dz;

            k1 = dFi(energy, mas[CurIndex], CurZ) * dz;
            k2 = dFi(energy, mas[CurIndex] + k1 / 2, CurZ + dz / 2) * dz;
            k3 = dFi(energy, mas[CurIndex] + k2 / 2, CurZ + dz / 2) * dz;
            k4 = dFi(energy, mas[CurIndex] + k3, CurZ + dz) * dz;

            mas[CurIndex + 1] = mas[CurIndex] + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
            fi[CurIndex] = mas[CurIndex];

            CurIndex++;
        }

        return mas[Controller.COUNT - 1];
    }

    double getEstat(double eMin, double eMax, double fiLeft, double fiRight, double statValue) {
        double midE;
        double midFi;
        int iter = 0;
        double err;

        do {
            midE = (eMax + eMin) / 2;
            midFi = getFi(midE);
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
        while (err > EPS);

        return midE;
    }

    double[] getR(double e) {

        double[] result = new double[COUNT];

        double zStart = -R;
        double zFinal = R;

        int CurIndex = COUNT / 2;
        double CurZ;
        double dz = (zFinal - zStart) / COUNT;
        double k1, k2, k3, k4;
        result[CurIndex] = 1;

        for (int i = CurIndex + 1; i < COUNT - 1; i++) {
            CurZ = zStart + i * dz;
            getFi(e);
            k1 = dz * dr(result[CurIndex], calculateU(coefficient, CurZ), e, fi[CurIndex]);
            k2 = dz * dr(result[CurIndex] + k1 / 2, calculateU(coefficient, CurZ + dz / 2), e, fi[CurIndex]);
            k3 = dz * dr(result[CurIndex] + k2 / 2, calculateU(coefficient, CurZ + dz / 2), e, fi[CurIndex]);
            k4 = dz * dr(result[CurIndex] + k3, calculateU(coefficient, CurZ + dz), e, fi[CurIndex]);

            result[CurIndex + 1] = result[CurIndex] + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
            CurIndex++;
        }

        dz = -(zFinal - zStart) / Controller.COUNT;
        CurIndex = Controller.COUNT / 2 - 1;
        result[CurIndex] = 1;

        for (int i = CurIndex - 1; i >= 0; i--) {
            CurZ = zStart - i * dz;
            getFi(e);
            k1 = dr(result[CurIndex], calculateU(coefficient, CurZ), e, fi[CurIndex]);
            k2 = dr(result[CurIndex] + k1 * dz / 2, calculateU(coefficient, CurZ + dz / 2), e, fi[CurIndex]);
            k3 = dr(result[CurIndex] + k2 * dz / 2, calculateU(coefficient, CurZ + dz / 2), e, fi[CurIndex]);
            k4 = dr(result[CurIndex] + k3 * dz, calculateU(coefficient, CurZ + dz), e, fi[CurIndex]);

            CurIndex--;
            result[CurIndex] = result[CurIndex + 1] + (k1 + 2 * k2 + 2 * k3 + k4) * dz / 6;

        }
        return result;
    }

    double dFi(double energy, double fi, double z) {
        return (calculateU(coefficient, z) - energy) * Math.pow(Math.cos(fi), 2) - Math.pow(Math.sin(fi), 2);
    }

    double dr(double r, double U, double E, double fi) {
        return r * (U - E + 1) * Math.cos(fi) * Math.sin(fi);
    }

    public void onDrawButton() {

        coefficient = Double.parseDouble(kText.getText());
        a = Double.parseDouble(pitWidthText.getText());
        U0 = Double.parseDouble(energyText.getText());
        sigma = Double.parseDouble(sigmaText.getText());

        fi = new double[COUNT];
        r = new double[COUNT];
        lastFi = new double[COUNT];

        double eMin = 0;
        double eMax = coefficient * a;

        double currentFi = -Math.PI / 2;

        List<Double> energyLevels = new ArrayList<>();

        lastFi[0] = getFi(eMin);
        for (int i = 1; i < COUNT; i++) {
            double currentE = eMin + ((eMax - eMin) / COUNT) * i;
            lastFi[i] = getFi(currentE);
            if (lastFi[i] < currentFi && lastFi[i - 1] > currentFi) {
                double energy = getEstat(currentE - (eMax - eMin) / COUNT, currentE, lastFi[i - 1], lastFi[i], currentFi);
                energyLevels.add(energy);
                currentFi = -(2 * energyLevels.size() + 1) * Math.PI / 2;
            }
        }

        // Draw phase chart
        XYChart.Series<Number, Number> series1 = new XYChart.Series<>();
        for (int i = 0; i < COUNT; i++) {
            series1.getData().add(new XYChart.Data<>(i, lastFi[i]));
        }
        phaseEnergyChart.getData().clear();
        phaseEnergyChart.getData().add(series1);

        // Draw wave function chart
        waveFunctionChart.getData().clear();
        for (Double energyLevel : energyLevels) {
            XYChart.Series<Number, Number> waveFunction = new XYChart.Series<>();
            r = getR(energyLevel);
            getFi(energyLevel);
            for (int j = 0; j < COUNT; j++) {
                double value = r[j] * Math.cos(fi[j]);
                waveFunction.getData().add(new XYChart.Data<>(j, value));
            }
            waveFunctionChart.getData().add(waveFunction);
        }
    }
}
