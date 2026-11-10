package scratch.anne.sra.gui;

import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Path2D;
import java.util.List;
import java.util.Map;
import java.util.ArrayList;
import java.util.Arrays;

import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.imr.AttenRelRef;

import scratch.anne.sra.domain.asset.fragility.FragilityAsset;
import scratch.anne.sra.domain.asset.fragility.PBRSurvivalAsset;
import scratch.anne.sra.domain.hazard.AssetHazardCurve;
import scratch.anne.sra.domain.hazard.HazardParameters;
import scratch.anne.sra.util.enums.HazardMetric;
import scratch.anne.sra.domain.portfolio.Portfolio;
import scratch.anne.sra.domain.portfolio.portfolioWrappers.PBRSurvivalPortfolioReporter;
import scratch.anne.sra.domain.portfolio.portfolioWrappers.RiskConvolutionPortfolio;
import scratch.anne.sra.domain.structuralResponse.ResponseModelLibrary;
import scratch.anne.sra.domain.structuralResponse.SimpleImResponse;
import scratch.anne.sra.domain.structuralResponse.SimpleImResponseLibrary;
import scratch.anne.sra.domain.structuralResponse.fragilities.FragilityBuilder;
import scratch.anne.sra.domain.structuralResponse.fragilities.FragilityModel;
import scratch.anne.sra.domain.structuralResponse.fragilities.SimpleImFragilityLibraryPreparer;
import scratch.anne.sra.domain.structuralResponse.fragilities.SimpleImFragilityLibraryPreparer.PrepareSpec;
import scratch.anne.sra.engine.convolution.ConvolutionResult;
import scratch.anne.sra.engine.portfolioWorkflow.PortfolioRiskConvolutionCalculator;
import scratch.anne.sra.util.Metadata;
import scratch.anne.sra.util.StringUtil;
import scratch.anne.sra.util.enums.IMT;
import scratch.anne.sra.util.enums.LimitState;


//TODO clean up PBR gui for updated workflow
//TODO make a gui based on the standard OpenSHA template

public class PBRAppUI extends JFrame {
	
    // =========================
    // UI STATE (FIELDS)
    // =========================

	// --- Site ---
	private JTextField latField, lonField, vs30Field;

	// --- Hazard ---
	private JComboBox<ERFOption> erfCombo;
	private JComboBox<GMMOption> gmmCombo;
	
	// --- Fragility ---
	private JTextField ageField, medianField, betaField;
	
	private JComboBox<IMTOption> imtCombo;

	// --- Probability ---
	private JTextField probField;
	
	// --- Output ---
	private PlotPanel plotPanel;
	private JTextArea resultsText;
	
	private JProgressBar progressBar;
	private JLabel progressLabel;

	
	private static final ERFOption[] ALLOWED_ERFS = {
		new ERFOption("org.opensha.sha.earthquake.rupForecastImpl.nshm23.erf.NSHM23_BranchAveragedERF"),
	    new ERFOption("org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final"),
	    new ERFOption("org.opensha.sha.earthquake.rupForecastImpl.MeanUCERF3"),
	};
	

	private static final GMMOption[] ALLOWED_GMMS = {
			new GMMOption(AttenRelRef.USGS_NSHM23_ACTIVE),
		    new GMMOption(AttenRelRef.ASK_2014),
		    new GMMOption(AttenRelRef.BSSA_2014),
		    new GMMOption(AttenRelRef.CB_2014)
		};
	
	private static final IMTOption[] ALLOWED_IMTS = {
		    new IMTOption(IMT.PGA, null),
		    new IMTOption(IMT.PGV, null),
		    new IMTOption(IMT.SA, 0.3),
		    new IMTOption(IMT.SA, 1.0),
		    new IMTOption(IMT.SA, 2.0)
		};


    public PBRAppUI() {
        super("Precarious Balanced Rock: Hazard Adjustment Calculator");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setSize(1000, 600);
        setLocationRelativeTo(null);

        setLayout(new BorderLayout());
        add(buildMainPanel(), BorderLayout.CENTER);
    }

    private JPanel buildMainPanel() {
        JPanel panel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc = defaultGbc();

        // LEFT: Plot
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.weightx = 0.6;
        gbc.weighty = 0.7;
        panel.add(createPlotPanel(), gbc);

        // LEFT: Results
        gbc.gridy = 1;
        gbc.weighty = 0.3;
        panel.add(createResultsPanel(), gbc);

        // RIGHT: Controls (scrollable)
        gbc.gridx = 1;
        gbc.gridy = 0;
        gbc.gridheight = 2;
        gbc.weightx = 0.4;
        gbc.weighty = 1.0;
        panel.add(buildControlPanel(), gbc);

        return panel;
    }

    private JScrollPane buildControlPanel() {
    	JPanel container = new JPanel();
    	container.setLayout(new BoxLayout(container, BoxLayout.Y_AXIS));

        JPanel formPanel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc = defaultGbc();

        gbc.weightx = 1;
        gbc.fill = GridBagConstraints.HORIZONTAL;

        gbc.gridy = 0;
        formPanel.add(createSitePanel(), gbc);

        gbc.gridy++;
        formPanel.add(createHazardPanel(), gbc);

        gbc.gridy++;
        formPanel.add(createFragilityPanel(), gbc);

        gbc.gridy++;
        formPanel.add(createProbabilityPanel(), gbc);

        // push everything up
        gbc.gridy++;
        gbc.weighty = 1;
        formPanel.add(Box.createVerticalGlue(), gbc);
        
        // Progress panel (above buttons)
        JPanel progressPanel = new JPanel(new BorderLayout());
        progressLabel = new JLabel("Idle");
        progressBar = new JProgressBar();
        progressBar.setIndeterminate(false);

        progressPanel.add(progressLabel, BorderLayout.NORTH);
        progressPanel.add(progressBar, BorderLayout.CENTER);
        progressPanel.setMaximumSize(
                new Dimension(Integer.MAX_VALUE, 45)
        );
        progressBar.setPreferredSize(new Dimension(200, 18));
        progressLabel.setBorder(
                BorderFactory.createEmptyBorder(2, 5, 2, 5)
        );

        // Bottom button panel
        JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.RIGHT, 5, 5));
        JButton calc = new JButton("Calculate");
        calc.addActionListener(e -> runCalculation());
        buttonPanel.add(calc);
        
        JButton clear = new JButton("Clear Results");
        clear.addActionListener(e -> {
            resultsText.setText("");
            plotPanel.clear();
        });
        buttonPanel.add(clear);

        container.add(formPanel);
        container.add(Box.createVerticalStrut(5));
        container.add(progressPanel);
        container.add(Box.createVerticalStrut(5));
        container.add(buttonPanel);

        JScrollPane scroll = new JScrollPane(container);
        scroll.setBorder(null);
        scroll.getVerticalScrollBar().setUnitIncrement(10);

        return scroll;
    }

    // ---------- Panels ----------

    private JPanel createPlotPanel() {
        plotPanel = new PlotPanel();
        plotPanel.setBorder(new TitledBorder(new EtchedBorder(), "Results Plot"));
        return plotPanel;
    }

    private JPanel createResultsPanel() {
        JPanel panel = new JPanel(new BorderLayout());
        panel.setBorder(new TitledBorder(new EtchedBorder(), "Results"));

        resultsText = new JTextArea();
        resultsText.setEditable(false);
        resultsText.setLineWrap(true);
        resultsText.setWrapStyleWord(true);

        JScrollPane scroll = new JScrollPane(resultsText);

        panel.add(scroll, BorderLayout.CENTER);
        return panel;
    }

    private JPanel createSitePanel() {
        JPanel panel = createSection("Site");
        
        latField = new JTextField("34.18", 10);
        lonField = new JTextField("-117.31", 10);
        vs30Field = new JTextField("760", 10);

        panel.add(new JLabel("Lat:"), labelGbc(0));
        panel.add(latField, fieldGbc(0));

        panel.add(new JLabel("Lon:"), labelGbc(1));
        panel.add(lonField, fieldGbc(1));

        panel.add(new JLabel("Vs30:"), labelGbc(2));
        panel.add(vs30Field, fieldGbc(2));

        return panel;
    }

    private JPanel createHazardPanel() {
        JPanel panel = createSection("Hazard");
        
        // ERF dropdown
        erfCombo = new JComboBox<>(ALLOWED_ERFS);

        // GMM dropdown
        gmmCombo = new JComboBox<>(ALLOWED_GMMS);

        panel.add(new JLabel("ERF:"), labelGbc(0));
        panel.add(erfCombo, fieldGbc(0));

        panel.add(new JLabel("GMM:"), labelGbc(1));
        panel.add(gmmCombo, fieldGbc(1));

        return panel;
    }

    private JPanel createFragilityPanel() {
        JPanel panel = createSection("Fragility");
        
        ageField = new JTextField("16", 10);
        medianField = new JTextField("1", 10);
        betaField = new JTextField("0.5", 10);
        
        imtCombo = new JComboBox<>(ALLOWED_IMTS);

        imtCombo.addActionListener(e -> {
            IMTOption sel = (IMTOption) imtCombo.getSelectedItem();
        });
        

        panel.add(new JLabel("Age (ka):"), labelGbc(0));
        panel.add(ageField, fieldGbc(0));

        panel.add(new JLabel("Median:"), labelGbc(1));
        panel.add(medianField, fieldGbc(1));

        panel.add(new JLabel("LogStdDev:"), labelGbc(2));
        panel.add(betaField, fieldGbc(2));

        panel.add(new JLabel("IMT:"), labelGbc(3));
        panel.add(imtCombo, fieldGbc(3));

        return panel;
    }

    private JPanel createProbabilityPanel() {
        JPanel panel = createSection("Probability of Survival");
        
        probField = new JTextField("0.05",10);

        panel.add(new JLabel("Target Probability:"), labelGbc(0));
        panel.add(probField, fieldGbc(0));

        return panel;
    }

    private JPanel createSection(String title) {
        JPanel panel = new JPanel(new GridBagLayout());
        panel.setBorder(new TitledBorder(new EtchedBorder(), title));
        panel.setOpaque(false);
        return panel;
    }

    // ---------- Layout Helpers ----------

    private GridBagConstraints defaultGbc() {
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(3, 3, 3, 3); // tighter OpenSHA-like spacing
        gbc.fill = GridBagConstraints.BOTH;
        return gbc;
    }

    private GridBagConstraints labelGbc(int y) {
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = y;
        gbc.anchor = GridBagConstraints.EAST;
        gbc.insets = new Insets(2, 2, 2, 4);
        return gbc;
    }

    private GridBagConstraints fieldGbc(int y) {
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.gridx = 1;
        gbc.gridy = y;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.weightx = 1;
        gbc.insets = new Insets(2, 2, 2, 2);
        return gbc;
    }

    // ---------- Main ----------

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            new PBRAppUI().setVisible(true);
        });
    }
    
    
    // --------- collect inputs and run ----------------
    private void runCalculation() {

        plotPanel.clear();
        
        setStage(Stage.IDLE);
        setStage(Stage.PARSING_INPUTS);

        SwingWorker<Void, Void> worker = new SwingWorker<>() {

            // store computed results as fields of the anonymous worker
            AssetHazardCurve hazard;
            FragilityModel baseFragility;
            SimpleImResponse fragility; 
            PBRSurvivalAsset asset;
            double hazardFactor;
            ConvolutionResult riskResult;

            @Override
            protected Void doInBackground() throws Exception {

                int assetID = 0;
                String fragilityID = "GUI_Fragility";
                
                setStage(Stage.PARSING_INPUTS);

                double lat = Double.parseDouble(latField.getText());
                double lon = Double.parseDouble(lonField.getText());
                double vs30 = Double.parseDouble(vs30Field.getText());

                double age = Double.parseDouble(ageField.getText());
                double median = Double.parseDouble(medianField.getText());
                double beta = Double.parseDouble(betaField.getText());
                double targetP = Double.parseDouble(probField.getText());

                ERFOption selectedERF = (ERFOption) erfCombo.getSelectedItem();
                GMMOption selectedGMM = (GMMOption) gmmCombo.getSelectedItem();
                IMTOption imtOption = (IMTOption) imtCombo.getSelectedItem();

                String erfClass = selectedERF.getClassName();
                AttenRelRef gmm = selectedGMM.getRef();
                String imtString = imtOption.getImtString();
                
                setStage(Stage.BUILDING_MODELS);

                baseFragility =
                        FragilityBuilder.component(fragilityID)
                                .imtString(imtString)
                                .lognormal(LimitState.TOPPLE, median, beta)
                                .build();

                ResponseModelLibrary<FragilityModel> baseFragilityLib =
                        ResponseModelLibrary.of(List.of(baseFragility), new Metadata.Builder().build());

                SimpleImResponseLibrary fragilityLib =
                        SimpleImFragilityLibraryPreparer.prepare(
                                baseFragilityLib,
                                PrepareSpec.builder().build()
                        );
                
                setStage(Stage.BUILDING_PORTFOLIO);

                List<FragilityAsset> assets = List.of(new FragilityAsset(
                        assetID, lat, lon, vs30, fragilityID,
                        Map.of("age", String.valueOf(age))
                ));

                Portfolio<FragilityAsset> basePortfolio =
                        new Portfolio<>(assets, List.of(), new Metadata.Builder().build());

                RiskConvolutionPortfolio riskPortfolio =
                        new RiskConvolutionPortfolio(basePortfolio, fragilityLib);
                
                HazardParameters hazardParameters = new HazardParameters(
                        erfClass,
                        1.0,
                        HazardMetric.PROBABILITY_EXCEEDANCE,
                		gmm.name()
                );

                // CAN NO LONGER SHOW THE ERF STEP OUTSIDE OF THE CALCULATOR
                //    but the .ASSESSING_RISK is [relatively] short so it just gets lumped with the BUILDING_ERF
                setStage(Stage.BUILDING_ERF);
//                
//                AbstractERF erf = buildERF(erfClass);
//                erf.getTimeSpan().setDuration(1.0);
//                erf.updateForecast();
//                
//                setStage(Stage.ASSESSING_RISK);
//
//                PortfolioRiskConvolutionCalculator calculator =
//                        new PortfolioRiskConvolutionCalculator(
//                                riskPortfolio,
//                                fragilityLib,
//                                gmm,
//                                erf
//                        );
                
                PortfolioRiskConvolutionCalculator calculator =
                        new PortfolioRiskConvolutionCalculator(
                                riskPortfolio,
                                fragilityLib,
                                hazardParameters
                        );

                riskPortfolio = calculator.computeRisk();

                setStage(Stage.POST_PROCESSING);
                
                PBRSurvivalPortfolioReporter reporter =
                        new PBRSurvivalPortfolioReporter(riskPortfolio);

                hazard = riskPortfolio.getHazardForAsset(assetID);
                riskResult = riskPortfolio.getAssetByID(assetID).getRiskConvolutionResult();
                
                fragility = fragilityLib.getModelByName(fragilityID);
                
                asset = reporter.getAssetByID(assetID);
                hazardFactor = asset.getHazardAdjustment(targetP);
                
                PlotData plotData = new PlotData(hazard, fragility, hazardFactor, riskResult);
                plotPanel.setPlotData(plotData);

                return null;
            }

            @Override
            protected void done() {

                // get() re-throws whatever doInBackground() threw, if
                // anything -- without this, SwingWorker silently swallows
                // that exception and done() proceeds anyway with
                // whatever fields doInBackground() managed to set
                // before failing (possibly none), leading to a confusing
                // downstream NPE instead of the real cause.
                try {
                    get();
                } catch (Exception e) {
                    setStage(Stage.IDLE);
                    Throwable cause = e.getCause() != null ? e.getCause() : e;
                    logResult(" ERROR: " + cause);
                    cause.printStackTrace();
                    return;
                }

            	setStage(Stage.COMPLETE);

                // ---------------- ORIGINAL LOG BLOCK (RESTORED) ----------------
                logResult(" New Calculation");
                logResult("================================");
                logResult(" Hazard");
                logResult("    ERF = " + erfCombo.getSelectedItem());
                logResult("    GMM = " + gmmCombo.getSelectedItem());
                logResult(String.format(
                        "    Lat, Lon = (%.5f, %.5f), Vs30 = %.0f",
                        Double.parseDouble(latField.getText()),
                        Double.parseDouble(lonField.getText()),
                        Double.parseDouble(vs30Field.getText())
                ));
                logResult("---------------------------------------------------------");
                logResult(" Fragility");
                logResult("    Age = " + ageField.getText() + "ka");
                logResult("    IMT = " + imtCombo.getSelectedItem());
                logResult("    Median = " + medianField.getText() + " [g]");
                logResult("    LogStdDev = " + betaField.getText());

                // ---------------- RESULTS ----------------
                logResult("---------------------------------------------------------");
                logResult(" Results");

                logResult(String.format(
                        "    Annual probability of failure = %.2e",
                        asset.getAnnualProbabilityOfFailure()
                ));

                logResult(String.format(
                        "    Probability of survival over %.0f ka = %.2e",
                        asset.getAge() / 1000,
                        asset.getProbabilityOfSurvival()
                ));

                logResult(String.format(
                        "    Hazard adjustment factor = %.2e",
                        hazardFactor
                ));
                
                logResult("");
                logResult("");
            }
        };

        worker.execute();
    }
    
    private AbstractERF buildERF(String className) throws Exception {

        Class<?> clazz = Class.forName(className);

        if (!AbstractERF.class.isAssignableFrom(clazz))
            throw new IllegalArgumentException(
                className + " is not an AbstractERF");

        return (AbstractERF)
                clazz.getDeclaredConstructor().newInstance();
    }
    
    private static class ERFOption {

        private final String className;
        private final String displayName;

        public ERFOption(String className) {
            this.className = className;
            this.displayName = className.substring(className.lastIndexOf('.') + 1);
        }

        public String getClassName() {
            return className;
        }

        @Override
        public String toString() {
            return displayName;
        }
    }
    
    private static class GMMOption {

        private final AttenRelRef ref;

        public GMMOption(AttenRelRef ref) {
            this.ref = ref;
        }

        public AttenRelRef getRef() {
            return ref;
        }

        @Override
        public String toString() {
            return ref.toString();
        }
    }
    
    private static class IMTOption {

        private final IMT imt;     // e.g. "SA", "PGA"
        private final Double period;  // null unless SA

        public IMTOption(IMT imt, Double period) {
            this.imt = imt;
            this.period = period;
        }

//        public Double getPeriod() {
//            return period;
//        }
        
        public String getImtString() {
        	return StringUtil.imtToString(imt, period);
        }
        
        @Override
        public String toString() {
            return StringUtil.imtToString(imt, period);
        }
    }
    
    
    public class FragilityAssetInput {
        public String assetID;
        public double lat;
        public double lon;
        public double vs30;
        public String responseModel;
        public double age;
    }

    
    private void logResult(String text) {
        resultsText.append(text + "\n");
        resultsText.setCaretPosition(resultsText.getDocument().getLength());
    }
    
    
    public class PlotData {

    	// hazard
        public final double[] imls;
        public final double[] poeOriginal;
        public final double hazardFactor;
        public final double[] poeAdjusted;
        
        // fragility
        public final double[] fragilityProbability;
        
        // disagg
        public final double[] disaggCumulative;
        public final double[] disaggBounds;

        public PlotData(AssetHazardCurve hazard,
        				SimpleImResponse fragility,
                        double hazFactor,
                        ConvolutionResult riskResult,
                        double[] disaggMarkerValues) {
        	
            this.imls = hazard.imls;
            this.poeOriginal = hazard.hazard;
            this.hazardFactor = hazFactor;

            this.poeAdjusted = new double[hazard.imls.length];
	            for (int i = 0; i < poeAdjusted.length; i++) {
	            	poeAdjusted[i] = hazard.hazard[i] * hazardFactor;
	            }

            this.fragilityProbability = fragility.getRespEdges();
            
            this.disaggCumulative = riskResult.getDisaggCumulative();
            this.disaggBounds = imlsAtMarkers(
            		imls,
            		disaggCumulative,
                    disaggMarkerValues);
        }
        
        public PlotData(AssetHazardCurve hazard,
				SimpleImResponse fragility,
                double hazFactor,
                ConvolutionResult riskResult
                ) {
            this(
                    hazard,
                    fragility,
                    hazFactor,
                    riskResult,
                    new double[] {0.25, 0.75}
                );
        }
        
        public static double imlAtProbability(
                double[] imls,
                double[] cdf,
                double targetProb) {

            int n = imls.length;

            int idx = Arrays.binarySearch(cdf, targetProb);

            if (idx >= 0) return imls[idx];

            int insert = -idx - 1;

            if (insert == 0) return imls[0];
            if (insert >= n) return imls[n - 1];

            double p0 = cdf[insert - 1];
            double p1 = cdf[insert];

            double x0 = imls[insert - 1];
            double x1 = imls[insert];

            double t = (targetProb - p0) / (p1 - p0);

            return x0 + t * (x1 - x0);
            

        }
        
        public static double[] imlsAtMarkers(
                double[] imls,
                double[] disaggCumulative,
                double[] markerValues) {

            double[] result = new double[markerValues.length];

            for (int m = 0; m < markerValues.length; m++) {

                double target = markerValues[m];

                result[m] = imlAtProbability(
                        imls,
                        disaggCumulative,
                        target);
            }

            return result;
        }  
    
    }
    
    
    private static class PlotPanel extends JPanel {

        /* ============================================================
         *                        AXIS MODEL
         * ============================================================ */

        public enum YAxis {
            LEFT_LOG,      // hazard curves (log Y)
            RIGHT_LINEAR   // fragility curves (linear Y)
        }
        
        private double[] verticalLines = new double[0];

        private static class Curve {
            double[] x;
            double[] y;
            String label;
            int index;
            YAxis yAxis;

            Color color; // <-- explicit override (null = use dynamic)

            Curve(double[] x,
                  double[] y,
                  String label,
                  int index,
                  YAxis yAxis,
                  Color color) {

                this.x = x;
                this.y = y;
                this.label = label;
                this.index = index;
                this.yAxis = yAxis;
                this.color = color;
            }
        }

        private final List<Curve> curves = new ArrayList<>();
        private PlotData data;

        /* ============================================================
         *                        STYLE
         * ============================================================ */

        private Color resolveColor(Curve c) {

            // explicit override wins
            if (c.color != null)
                return c.color;

            // otherwise dynamic palette
            return Color.getHSBColor(c.index / 10f, 0.8f, 0.9f);
        }

        private final int leftMargin   = 80;
        private final int rightMargin  = 60;
        private final int topMargin    = 40;
        private final int bottomMargin = 50;

        public PlotPanel() {
        }

        /* ============================================================
         *                        PUBLIC API
         * ============================================================ */

        public void setPlotData(PlotData data) {
            this.data = data;
            this.verticalLines = data.disaggBounds;

            curves.clear();
            
            
            // fragility curve (RIGHT axis, linear)
            curves.add(new Curve(
                    data.imls,
                    data.fragilityProbability,
                    "Fragility",
                    2,
                    YAxis.RIGHT_LINEAR,
                    Color.BLACK   // fixed semantic color
            ));

            // hazard curves (LEFT axis, log)
            curves.add(new Curve(
                    data.imls,
                    data.poeOriginal,
                    "Computed Hazard",
                    0,
                    YAxis.LEFT_LOG,
                    null   // dynamic
            ));

            curves.add(new Curve(
                    data.imls,
                    data.poeAdjusted,
                    String.format("Adjusted Hazard (×%.2g)", data.hazardFactor),
                    1,
                    YAxis.LEFT_LOG,
                    null   // dynamic
            ));
            

            repaint();
        }

        public void clear() {
            curves.clear();
            data = null;
            repaint();
        }

        /* ============================================================
         *                        PAINT
         * ============================================================ */

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);

            if (data == null || curves.isEmpty())
                return;

            Graphics2D g2 = (Graphics2D) g.create();
            g2.setRenderingHint(
                    RenderingHints.KEY_ANTIALIASING,
                    RenderingHints.VALUE_ANTIALIAS_ON
            );

            int w = getWidth();
            int h = getHeight();

            int L = leftMargin;
            int R = w - rightMargin;
            int T = topMargin;
            int B = h - bottomMargin;

            Limits lim = computeLimits();

            drawAxes(g2, L, R, T, B);
            drawTicks(g2, L, R, T, B, lim);

            /* ---------------- curves ---------------- */

            for (Curve c : curves) {
            	g2.setColor(resolveColor(c));
                drawCurve(g2, c, L, R, T, B, lim);
            }
            
            
            /* ---------------- vertical lines ---------------- */
            for (double x : this.verticalLines) {
                drawVertical(g2, x, L, R, T, B, lim);
            }

            /* ---------------- legend ---------------- */

            drawLegend(g2, L, R, T, B);

            /* ---------------- labels ---------------- */

            drawLeftLabel(g2, "Annual Probability of Exceedance", L, T, B);
            drawRightLabel(g2, "P (Collapse | IML)", R, T, B, w);
            drawXLabel(g2,
                    "Intensity Measure Level (g)",
                    L, R, B);

            g2.dispose();
        }

        /* ============================================================
         *                        CURVE DRAW
         * ============================================================ */

        private void drawCurve(Graphics2D g,
                               Curve c,
                               int L, int R, int T, int B,
                               Limits lim) {

			            Path2D path = new Path2D.Double();
			
			            for (int i = 0; i < c.x.length; i++) {
			
			                int px = mapLogX(c.x[i], L, R, lim);
			
			                int py;
			                if (c.yAxis == YAxis.LEFT_LOG) {
			                    py = mapLogY(c.y[i], T, B, lim);
			                } else {
			                    py = mapLinearY(c.y[i], T, B);
			                }
			
			                if (i == 0)
			                    path.moveTo(px, py);
			                else
			                    path.lineTo(px, py);
			            }
			
			            g.draw(path);
			        }
			        
			        
			        private void drawVertical(Graphics2D g,
			                double x,
			                int L, int R, int T, int B,
			                Limits lim) {
			
						if (data == null) return;
						if (!Double.isFinite(x) || x <= 0) return;
						
						double logX = Math.log10(x);
						
						double denom = (lim.logXmax - lim.logXmin);
						if (denom == 0) return;
						
						double t = (logX - lim.logXmin) / denom;
						int px = (int) (L + t * (R - L));
						
						Stroke oldStroke = g.getStroke();
						Color oldColor = g.getColor();
						AffineTransform oldTransform = g.getTransform();

						FontMetrics fm = g.getFontMetrics();

						// ---------------- vertical line ----------------
						g.setColor(Color.GRAY);

						g.setStroke(new BasicStroke(
						        1.5f,
						        BasicStroke.CAP_BUTT,
						        BasicStroke.JOIN_BEVEL,
						        0,
						        new float[]{6f, 6f},
						        0));

						g.drawLine(px, T, px, B);

						// ---------------- label ----------------
						String label = String.format("%.2g", x);
						int textW = fm.stringWidth(label);

						int xDraw = px + 3; // - textW / 2;
						int yDraw = B - fm.getDescent() - 3;

						g.drawString(label, xDraw, yDraw);

						// ---------------- restore ALL STATE ONCE ----------------
						g.setTransform(oldTransform);
						g.setStroke(oldStroke);
						g.setColor(oldColor);
		}

        /* ============================================================
         *                        AXES
         * ============================================================ */

        private void drawAxes(Graphics2D g, int L, int R, int T, int B) {

            g.setColor(Color.BLACK);

            // bottom x-axis
            g.drawLine(L, B, R, B);

            // left y-axis (hazard)
            g.drawLine(L, T, L, B);

            // right y-axis (fragility)
            g.drawLine(R, T, R, B);
        }

        /* ============================================================
         *                        TICKS
         * ============================================================ */

        private void drawTicks(Graphics2D g,
                               int L, int R, int T, int B,
                               Limits lim) {

            FontMetrics fm = g.getFontMetrics();

            /* -------- X (log) -------- */

            for (int p = (int) Math.floor(lim.logXmin);
                    p <= (int) Math.ceil(lim.logXmax);
                    p++) {

                double base = Math.pow(10, p);

                // minor + major ticks per decade
                for (int m = 1; m <= 9; m++) {

                    double val = m * base;

                    if (val < lim.xmin || val > lim.xmax)
                        continue;

                    int x = mapLogX(val, L, R, lim);

                    boolean major = (m == 1);

                    // tick length
                    int tickLen = major ? 7 : 4;

                    g.drawLine(x, B, x, B + tickLen);

                    // label ONLY major ticks (10^p)
                    if (major) {
                        String label = "10^" + p;
                        int w = fm.stringWidth(label);
                        g.drawString(label, x - w / 2, B + 20);
                    }
                }
            }

            /* -------- LEFT Y (log hazard) -------- */

            for (int p = (int) lim.logYmin; p <= (int) lim.logYmax; p++) {

                double val = Math.pow(10, p);
                int y = mapLogY(val, T, B, lim);

                g.drawLine(L - 5, y, L, y);

                String label = "10^" + p;
                int w = fm.stringWidth(label);

                g.drawString(label, L - 10 - w, y + fm.getAscent() / 2);
            }

            /* -------- RIGHT Y (linear fragility) -------- */

            int steps = 5; // 0.0 → 1.0

            for (int i = 0; i <= steps; i++) {

                double val = i / (double) steps;
                int y = mapLinearY(val, T, B);

                // tick mark on right axis
                g.drawLine(R, y, R + 5, y);

                String label = String.format("%.1f", val);
                int w = fm.stringWidth(label);

                // draw label to the right of axis
                g.drawString(label, R + 8, y + fm.getAscent() / 2 - 2);
            }
        }

        /* ============================================================
         *                        LEGEND
         * ============================================================ */

        private void drawLegend(Graphics2D g, int L, int R, int T, int B) {

            FontMetrics fm = g.getFontMetrics();
            int lineH = fm.getHeight();

            int midY = (T + B) / 2;

            // center the whole legend block vertically
            int startY = midY - (curves.size() * lineH) / 2;

            for (int i = 0; i < curves.size(); i++) {

                Curve c = curves.get(i);

                g.setColor(resolveColor(c));

                int y = startY + i * lineH;

                g.drawString(
                        c.label,
                        L + 10,
                        y
                );
            }
        }

        /* ============================================================
         *                        LABELS
         * ============================================================ */

        private void drawLeftLabel(Graphics2D g,
			                String text,
			                int L,
			                int T,
			                int B) {
			
			AffineTransform old = g.getTransform();
			FontMetrics fm = g.getFontMetrics();
			
			g.setColor(Color.BLACK);
			
			g.rotate(-Math.PI / 2);
			
			int textWidth = fm.stringWidth(text);
			int centerY = T + (B - T) / 2;
			
			// IMPORTANT: push outside plot area, not onto axis
			int x = -(centerY + textWidth / 2);
			int y = L - 50;   // <-- increased gutter beyond tick labels
			
			g.drawString(text, x, y);
			
			g.setTransform(old);
		}

        private void drawRightLabel(Graphics2D g,
                String text,
                int R,
                int T,
                int B,
                int w) {

			AffineTransform old = g.getTransform();
			FontMetrics fm = g.getFontMetrics();
			
			g.setColor(Color.BLACK);
			
			g.rotate(-Math.PI / 2);
			
			int textWidth = fm.stringWidth(text);
			int centerY = T + (B - T) / 2;
			
			int x = -(centerY + textWidth / 2);
			int y = R + 40;  
			
			g.drawString(text, x, y);
			
			g.setTransform(old);
		}
        
        private void drawXLabel(Graphics2D g,
			                String text,
			                int L, int R, int B) {
			
			FontMetrics fm = g.getFontMetrics();
			
			int textWidth = fm.stringWidth(text);
			
			int x = L + (R - L) / 2 - textWidth / 2;
			int y = B + fm.getAscent() + 25;  // padding below axis
			
			g.setColor(Color.BLACK);
			g.drawString(text, x, y);
			}

        /* ============================================================
         *                        MAPPERS
         * ============================================================ */

        private int mapLogX(double x, int L, int R, Limits lim) {
            double t = (Math.log10(x) - lim.logXmin) /
                    (lim.logXmax - lim.logXmin);
            return (int) (L + t * (R - L));
        }

        private int mapLogY(double y, int T, int B, Limits lim) {
            double t = (Math.log10(y) - lim.logYmin) /
                    (lim.logYmax - lim.logYmin);
            return (int) (B - t * (B - T));
        }

        private int mapLinearY(double y, int T, int B) {
            return (int) (B - y * (B - T));
        }

        /* ============================================================
         *                        LIMITS
         * ============================================================ */

        private Limits computeLimits() {

            Limits lim = new Limits();

            expand(lim, data.imls, data.poeOriginal);
            expand(lim, data.imls, data.poeAdjusted);

            lim.finalizeLogs();
            return lim;
        }

        private void expand(Limits lim, double[] x, double[] y) {

            if (x == null || y == null) return;

            for (int i = 0; i < x.length; i++) {

                lim.xmin = Math.min(lim.xmin, x[i]);
                lim.xmax = Math.max(lim.xmax, x[i]);
                lim.ymin = Math.min(lim.ymin, y[i]);
                lim.ymax = Math.max(lim.ymax, y[i]);
            }
        }

        /* ============================================================
         *                        LIMIT TYPE
         * ============================================================ */

        private static class Limits {

            double xmin = Double.POSITIVE_INFINITY;
            double xmax = Double.NEGATIVE_INFINITY;
            double ymin = Double.POSITIVE_INFINITY;
            double ymax = Double.NEGATIVE_INFINITY;

            double logXmin, logXmax;
            double logYmin, logYmax;

            void finalizeLogs() {
                logXmin = Math.log10(xmin);
                logXmax = Math.log10(xmax);
                logYmin = Math.log10(ymin);
                logYmax = Math.log10(ymax);
            }
        }
    }
    
    
    
    private void setStage(Stage stage) {
        progressLabel.setText(stage.label);

        switch (stage) {
            case IDLE -> {
                progressBar.setIndeterminate(false);
                progressBar.setValue(0);
            }
            case PARSING_INPUTS -> progressBar.setIndeterminate(true);
            case BUILDING_MODELS -> progressBar.setIndeterminate(true);
            case BUILDING_PORTFOLIO -> progressBar.setIndeterminate(true);
            case BUILDING_ERF -> progressBar.setIndeterminate(true);
            case ASSESSING_RISK -> progressBar.setIndeterminate(true);
            case POST_PROCESSING -> progressBar.setIndeterminate(true);
            case COMPLETE -> {
                progressBar.setIndeterminate(false);
                progressBar.setValue(100);
            }
        }
    }
    
    private enum Stage {
        IDLE("Idle"),
        PARSING_INPUTS("Reading inputs..."),
        BUILDING_MODELS("Building models..."),
        BUILDING_PORTFOLIO("Building portfolio..."),
        BUILDING_ERF("Building ERF..."),
        ASSESSING_RISK("Assessing risk..."),
        POST_PROCESSING("Post-processing..."),
        COMPLETE("Complete");

        final String label;

        Stage(String label) {
            this.label = label;
        }
    }

}
    

