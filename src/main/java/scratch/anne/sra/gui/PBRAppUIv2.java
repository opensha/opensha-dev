package scratch.anne.sra.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.lang.reflect.InvocationTargetException;
import java.util.List;
import java.util.Map;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.border.EtchedBorder;
import javax.swing.border.TitledBorder;

import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.gui.plot.GraphWidget;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.sha.earthquake.BaseERF;
import org.opensha.sha.earthquake.ERF_Ref;
import org.opensha.sha.gui.beans.ERF_GuiBean;
import org.opensha.sha.gui.beans.IMR_GuiBean;
import org.opensha.sha.gui.beans.IMR_GuiBeanAPI;
import org.opensha.sha.gui.beans.IMT_NewGuiBean;
import org.opensha.sha.gui.beans.Site_GuiBean;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;

import scratch.anne.sra.domain.asset.fragility.FragilityAsset;
import scratch.anne.sra.domain.asset.fragility.PBRSurvivalAsset;
import scratch.anne.sra.domain.hazard.AssetHazardCurve;
import scratch.anne.sra.domain.hazard.HazardParameters;
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
import scratch.anne.sra.util.enums.HazardMetric;
import scratch.anne.sra.util.enums.LimitState;

/**
 * PBR hazard-adjustment GUI, rebuilt on OpenSHA's standard
 * {@code *_GuiBean} architecture (see {@code HazardCurveApplication}
 * for the canonical example) in place of the hand-rolled combo boxes /
 * manual text parsing in the original {@code PBRAppUI}.
 *
 * This is a mockup, as the calculation pipeline is not yet set up to 
 * work this type of ERF or GMM creation. Instead, it is built inside
 * the calculator based on input strings and an erf.config file (limited
 * to sra package, not a standard format). Thus, it only takes the GMM
 * and ERFs name and builds them as their default version.
 * 
 * The portfolio assets assume a Vs30 value (consistent with Keith Porter's
 * portfolio csv files, so any GMM that doesn't use Vs30 is invalid,
 * since the Vs30 field goes away and cannot be read into an asset.
 */
public class PBRAppUIv2 extends JFrame implements IMR_GuiBeanAPI {

    private static final long serialVersionUID = 1L;

    // ------------------------------------------------------------
    // Standard OpenSHA hazard-setup beans
    // ------------------------------------------------------------

    private ERF_GuiBean erfGuiBean;
    private IMR_GuiBean imrGuiBean;
    private Site_GuiBean siteGuiBean;
    private IMT_NewGuiBean imtGuiBean;

    // ------------------------------------------------------------
    // PBR-specific inputs (ported from the original PBRAppUI --
    // no OpenSHA bean equivalent exists for these)
    // ------------------------------------------------------------

    private JTextField ageField, medianField, betaField;
    private JTextField probField;

    // ------------------------------------------------------------
    // Output (ported from the original PBRAppUI -- plot not yet wired)
    // ------------------------------------------------------------

    private JTextArea resultsText;
    private GraphWidget plotWidget;

    public PBRAppUIv2() {
        super("Precarious Balanced Rock: Hazard Adjustment Calculator (v2 mockup)");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setSize(1100, 700);
        setLocationRelativeTo(null);

        initBeans();

        setLayout(new BorderLayout());
        JSplitPane split = new JSplitPane(
                JSplitPane.HORIZONTAL_SPLIT,
                buildOutputPanel(),
                buildControlsScrollPane());
        split.setResizeWeight(0.55);
        add(split, BorderLayout.CENTER);
    }

    /**
     * Builds the 4 hazard-setup beans. Order matters: the IMT bean is
     * built from the IMR bean's currently-selected GMM, and the site
     * bean's parameters come from that same GMM's own
     * {@code getSiteParams()}.
     */
    private void initBeans() {
        try {
            // lists every registered ERF -- swap/add ERF_Ref values here
            // to change what's offered, no combo-box code required
            erfGuiBean = new ERF_GuiBean(
                    ERF_Ref.get(true, org.opensha.commons.util.ServerPrefUtils.SERVER_PREFS));
        } catch (InvocationTargetException e) {
            throw new RuntimeException("Failed to build ERF_GuiBean", e);
        }

        imrGuiBean = new IMR_GuiBean(this);

        ScalarIMR selectedIMR = imrGuiBean.getSelectedIMR_Instance();
        imtGuiBean = new IMT_NewGuiBean(selectedIMR);

        siteGuiBean = new Site_GuiBean();
        siteGuiBean.addSiteParams(selectedIMR.getSiteParams());
    }

    // ------------------------------------------------------------
    // IMR_GuiBeanAPI -- required hooks, called by imrGuiBean whenever
    // the selected GMM changes
    // ------------------------------------------------------------

    @Override
    public void updateIM() {
        // rebuild the IMT bean against the newly-selected GMM, so the
        // offered IMTs (and, for SA, the offered periods) always match
        // what that GMM actually supports
        ScalarIMR selectedIMR = imrGuiBean.getSelectedIMR_Instance();
        imtGuiBean.setIMR(selectedIMR);
    }

    @Override
    public void updateSiteParams() {
        // replace the site bean's parameter list with the newly-
        // selected GMM's own getSiteParams() (lat/lon are preserved --
        // see Site_GuiBean.replaceSiteParams)
        ScalarIMR selectedIMR = imrGuiBean.getSelectedIMR_Instance();
        siteGuiBean.replaceSiteParams(selectedIMR.getSiteParams());
    }

    // ------------------------------------------------------------
    // Layout
    // ------------------------------------------------------------

    private JPanel buildOutputPanel() {
        JPanel panel = new JPanel(new BorderLayout());

        plotWidget = new GraphWidget();
        plotWidget.setPreferredSize(new Dimension(600, 450));

        // GraphWidget's built-in "Log scale: X / Y" checkboxes are single,
        // shared controls -- they can't express "log for the hazard
        // subplot, linear for the fragility subplot" at all, and clicking
        // either one calls GraphWidget.setY_Log(...)/setX_Log(...), which
        // redraws EVERY subplot with that one shared flag (see
        // GraphWidget.drawGraph()), silently overwriting the correct
        // per-subplot axes updatePlot() sets up. Since there's no way for
        // this control to produce a correct state for this plot, hide it
        // rather than leave a control that can only ever be wrong.
        plotWidget.getButtonControlPanel().getCheckboxRow().setVisible(false);

        resultsText = new JTextArea();
        resultsText.setEditable(false);
        resultsText.setBorder(new TitledBorder(new EtchedBorder(), "Results"));
        JScrollPane resultsScroll = new JScrollPane(resultsText);
        resultsScroll.setPreferredSize(new Dimension(600, 200));

        JSplitPane plotResultsSplit = new JSplitPane(
                JSplitPane.VERTICAL_SPLIT, true, plotWidget, resultsScroll);
        plotResultsSplit.setResizeWeight(0.7);

        panel.add(plotResultsSplit, BorderLayout.CENTER);
        return panel;
    }

    /**
     * Stacks the 4 hazard-setup beans in resizable split panes (same
     * pattern as {@code HazardCurveApplication}'s {@code imrImtSplitPane}/
     * {@code imrImtSiteSplitPane}), then the PBR-specific fields and
     * the Calculate button below, all inside one scroll pane -- so a
     * bean that grows at runtime (the IMT bean's period row) or a
     * window that's too short to show everything at once both remain
     * reachable instead of being silently clipped.
     */
    private JScrollPane buildControlsScrollPane() {
        JPanel container = new JPanel();
        container.setLayout(new BoxLayout(container, BoxLayout.Y_AXIS));

        JSplitPane imrImtSplit = new JSplitPane(
                JSplitPane.VERTICAL_SPLIT, true, imrGuiBean, imtGuiBean);
        imrImtSplit.setResizeWeight(0.8);
        imrImtSplit.setPreferredSize(new Dimension(420, 320));

        JSplitPane imrImtSiteSplit = new JSplitPane(
                JSplitPane.VERTICAL_SPLIT, true, imrImtSplit, siteGuiBean);
        imrImtSiteSplit.setResizeWeight(0.75);
        imrImtSiteSplit.setPreferredSize(new Dimension(420, 560));

        JSplitPane erfSplit = new JSplitPane(
                JSplitPane.VERTICAL_SPLIT, true, imrImtSiteSplit, erfGuiBean);
        erfSplit.setResizeWeight(0.8);
        erfSplit.setPreferredSize(new Dimension(420, 700));

        // erfSplit (and the bean column generally) has a fixed preferred
        // width; when the controls column is wider than that (e.g. after
        // dragging the main split wider), BoxLayout.Y_AXIS centers each
        // child by default rather than left-aligning it, leaving a gap on
        // the LEFT of the fixed-width bean column that can't be closed by
        // resizing the plot box. Rather than fight BoxLayout's alignment,
        // wrap the beans in a BorderLayout (beans at EAST, disclaimer at
        // CENTER) so that gap is filled with something useful instead of
        // sitting empty.
        JPanel beansWithDisclaimer = new JPanel(new BorderLayout());
        beansWithDisclaimer.add(buildDisclaimerLabel(), BorderLayout.CENTER);
        beansWithDisclaimer.add(erfSplit, BorderLayout.EAST);

        container.add(beansWithDisclaimer);
        container.add(createFragilityPanel());
        container.add(createProbabilityPanel());

        JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.RIGHT, 5, 5));
        JButton calc = new JButton("Calculate");
        calc.addActionListener(e -> runCalculation());
        buttonPanel.add(calc);
        container.add(buttonPanel);

        JScrollPane scroll = new JScrollPane(container);
        scroll.getVerticalScrollBar().setUnitIncrement(12);
        scroll.setPreferredSize(new Dimension(460, 700));
        return scroll;
    }

    /**
     * Static disclaimer filling the fixed gap that {@code BoxLayout}
     * leaves to the left of the fixed-width bean column (see
     * {@link #buildControlsScrollPane()}) -- not a popup, since this is
     * a standing limitation of the current mockup, not a one-time
     * warning tied to a specific user action (contrast with the Vs30
     * check in {@link #runCalculation()}, which IS a popup because it's
     * about one specific attempted calculation).
     */
    private JTextArea buildDisclaimerLabel() {
        JTextArea disclaimer = new JTextArea(
                "While the IMR, Site, and ERF beans "
                + "are shown here, "
                + "not all of their features are yet included in "
                + "the calculation pipeline.\n"
                + "\n"
                + "Similarly, Vs30 is required "
                + "such that GMMs that don't use "
                + "Vs30 will not run."
        );
        disclaimer.setEditable(false);
        disclaimer.setFocusable(false);
        disclaimer.setOpaque(false);
        disclaimer.setLineWrap(true);
        disclaimer.setWrapStyleWord(true);
        disclaimer.setForeground(Color.DARK_GRAY);
        disclaimer.setFont(disclaimer.getFont().deriveFont(java.awt.Font.ITALIC));
        disclaimer.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
        return disclaimer;
    }

    private JPanel createFragilityPanel() {
        JPanel panel = new JPanel(new GridBagLayout());
        panel.setBorder(new TitledBorder(new EtchedBorder(), "Fragility"));
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(2, 4, 2, 4);
        gbc.anchor = GridBagConstraints.WEST;

        ageField = new JTextField("16", 10);
        medianField = new JTextField("1", 10);
        betaField = new JTextField("0.5", 10);

        gbc.gridx = 0; gbc.gridy = 0; panel.add(new JLabel("Age (ka):"), gbc);
        gbc.gridx = 1; panel.add(ageField, gbc);

        gbc.gridx = 0; gbc.gridy = 1; panel.add(new JLabel("Median:"), gbc);
        gbc.gridx = 1; panel.add(medianField, gbc);

        gbc.gridx = 0; gbc.gridy = 2; panel.add(new JLabel("LogStdDev:"), gbc);
        gbc.gridx = 1; panel.add(betaField, gbc);

        return panel;
    }

    private JPanel createProbabilityPanel() {
        JPanel panel = new JPanel(new GridBagLayout());
        panel.setBorder(new TitledBorder(new EtchedBorder(), "Probability of Survival"));
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(2, 4, 2, 4);
        gbc.anchor = GridBagConstraints.WEST;

        probField = new JTextField("0.05", 10);

        gbc.gridx = 0; gbc.gridy = 0; panel.add(new JLabel("Target Probability:"), gbc);
        gbc.gridx = 1; panel.add(probField, gbc);

        return panel;
    }

    // ------------------------------------------------------------
    // Calculation
    // ------------------------------------------------------------

    /**
     * Reads each bean's current selection on demand -- same
     * button-driven pattern as the original {@code PBRAppUI}, since
     * none of these beans push live change notifications to an outer
     * frame -- then runs the same fragility-build /
     * risk-convolution / hazard-adjustment pipeline the original did.
     *
     * <p>Unlike the original, the ERF class name and GMM
     * {@link AttenRelRef} aren't read from a custom combo-box option
     * type; they're derived from the beans' actual selected instances
     * ({@link BaseERF#getClass()}, and an {@link AttenRelRef} lookup
     * by matching {@link ScalarIMR#getName()} -- see
     * {@link #findAttenRelRef(ScalarIMR)}) -- {@code HazardParameters}/
     * {@code PortfolioRiskConvolutionCalculator} still only accept
     * string/enum identifiers and always rebuild their own ERF/GMM
     * from them, so those identifiers still have to be derived, but
     * now from whatever the user actually picked rather than from a
     * hardcoded option list.</p>
     * 
     */
    // TODO PBR gui's gmm and erf beans aren't fully ingested.
    private void runCalculation() {
        try {
            BaseERF erf = erfGuiBean.getSelectedERF();
            ScalarIMR imr = imrGuiBean.getSelectedIMR_Instance();
            Site site = siteGuiBean.getSite();
            String imtString = imtGuiBean.getSelectedIMT();

            double age = Double.parseDouble(ageField.getText());
            double median = Double.parseDouble(medianField.getText());
            double beta = Double.parseDouble(betaField.getText());
            double targetP = Double.parseDouble(probField.getText());

            String erfClass = erf.getClass().getName();
            AttenRelRef gmm = findAttenRelRef(imr);

            double lat = site.getLocation().getLatitude();
            double lon = site.getLocation().getLongitude();

            if (!site.containsParameter(Vs30_Param.NAME)) {
                javax.swing.JOptionPane.showMessageDialog(this,
                        "The asset class is requires a Vs30 parameter."
                        + "The selected GMM (" + imr.getName() + ") "
                        + "does not use a Vs30 site parameter."
                        + " Choose a different GMM.",
                        "Vs30 not available for this GMM",
                        javax.swing.JOptionPane.WARNING_MESSAGE);
                return;
            }
            double vs30 = (Double) site.getParameter(Vs30_Param.NAME).getValue();

            int assetID = 0;
            String fragilityID = "GUI_Fragility";

            FragilityModel baseFragility =
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

            PortfolioRiskConvolutionCalculator calculator =
                    new PortfolioRiskConvolutionCalculator(
                            riskPortfolio,
                            fragilityLib,
                            hazardParameters
                    );

            riskPortfolio = calculator.computeRisk();

            PBRSurvivalPortfolioReporter reporter =
                    new PBRSurvivalPortfolioReporter(riskPortfolio);

            PBRSurvivalAsset asset = reporter.getAssetByID(assetID);
            double hazardFactor = asset.getHazardAdjustment(targetP);

            AssetHazardCurve hazard = riskPortfolio.getHazardForAsset(assetID);
            SimpleImResponse fragilityResponse = fragilityLib.getModelByName(fragilityID);
            ConvolutionResult riskResult =
                    riskPortfolio.getAssetByID(assetID).getRiskConvolutionResult();

            updatePlot(hazard, fragilityResponse, hazardFactor, riskResult);

            logResult(" New Calculation");
            logResult("================================");
            logResult(" Hazard");
            logResult("    ERF = " + erfClass);
            logResult("    GMM = " + imr.getName());
            logResult(String.format("    Lat, Lon = (%.5f, %.5f), Vs30 = %.0f", lat, lon, vs30));
            logResult("---------------------------------------------------------");
            logResult(" Fragility");
            logResult("    Age = " + ageField.getText() + "ka");
            logResult("    IMT = " + imtString);
            logResult("    Median = " + medianField.getText() + " [g]");
            logResult("    LogStdDev = " + betaField.getText());
            logResult("---------------------------------------------------------");
            logResult(" Results");
            logResult(String.format("    Annual probability of failure = %.2e",
                    asset.getAnnualProbabilityOfFailure()));
            logResult(String.format("    Probability of survival over %.0f ka = %.2e",
                    asset.getAge() / 1000.0, asset.getProbabilityOfSurvival()));
            logResult(String.format("    Hazard adjustment factor = %.2e", hazardFactor));
            logResult("");

        } catch (InvocationTargetException e) {
            throw new RuntimeException("Failed to get selected ERF", e);
        }
    }

    /**
     * Finds the {@link AttenRelRef} whose {@code get().getName()}
     * matches {@code imr}'s own name -- {@link IMR_GuiBean} exposes
     * the selected GMM as a built {@link ScalarIMR} instance, not the
     * {@link AttenRelRef} identifier {@code HazardParameters} needs,
     * so this recovers it by name rather than requiring a second,
     * separately-maintained GMM list.
     */
    private AttenRelRef findAttenRelRef(ScalarIMR imr) {
        for (AttenRelRef ref : AttenRelRef.values()) {
            if (ref.getName().equals(imr.getName()))
                return ref;
        }
        throw new IllegalStateException(
                "No AttenRelRef found matching selected IMR name: " + imr.getName());
    }

    // ------------------------------------------------------------
    // Plot
    // ------------------------------------------------------------

    /**
     * Rebuilds {@link #plotWidget} as 2 stacked subplots sharing one
     * X (IML) axis -- hazard (log-Y) on top, fragility (linear-Y)
     * below -- via {@code GraphPanel}'s own {@code CombinedDomainXYPlot}
     * support, reached directly through {@link #plotWidget}'s
     * {@code getGraphPanel()} since {@link GraphWidget}'s own
     * convenience constructors/setters only take one shared log-axis
     * flag for every subplot, not a per-subplot flag as needed here.
     *
     * <p>The 2 disaggregation-bound vertical lines are added as
     * {@link ValueMarker}s on every subplot individually (not just
     * one) -- {@code CombinedDomainXYPlot} shares its domain axis
     * range across subplots, but a marker added to one subplot is
     * NOT automatically drawn on the others; each subplot's own
     * {@code addDomainMarker} must be called so the lines visually
     * span the full stacked figure.</p>
     */
    // TODO the stacked plot goes away when changing plot preferences
    private void updatePlot(AssetHazardCurve hazard, SimpleImResponse fragility,
                             double hazardFactor, ConvolutionResult riskResult) {

        DiscretizedFunc hazardOriginal = buildFunc(hazard.imls, hazard.hazard, "Computed Hazard");

        double[] poeAdjusted = new double[hazard.hazard.length];
        for (int i = 0; i < poeAdjusted.length; i++)
            poeAdjusted[i] = hazard.hazard[i] * hazardFactor;
        DiscretizedFunc hazardAdjusted = buildFunc(hazard.imls, poeAdjusted,
                String.format("Adjusted Hazard (\u00d7%.2g)", hazardFactor));

        DiscretizedFunc fragilityCurve = buildFunc(
                fragility.getImValues(), fragility.getRespEdges(), "Fragility");

        PlotSpec hazardSpec = new PlotSpec(
                List.of(hazardOriginal, hazardAdjusted),
                List.of(
                        new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE),
                        new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED)),
                "", "IML", "Probability of Exceedance");
        hazardSpec.setLegendVisible(true);

        PlotSpec fragilitySpec = new PlotSpec(
                List.of(fragilityCurve),
                List.of(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK)),
                "", "IML", "Fragility Probability");
        fragilitySpec.setLegendVisible(true);

        double[] disaggBounds = imlsAtMarkers(hazard.imls, riskResult.getDisaggCumulative(),
                new double[] {0.25, 0.75});

        // GraphWidget's own convenience methods (setPlotSpec/
        // setMultiplePlotSpecs/drawGraph) only support one shared
        // xLog/yLog flag applied to every subplot -- they can't express
        // "log for the hazard subplot, linear for the fragility subplot"
        // at all. So this calls GraphPanel's own lower-level
        // drawGraphPanel(...) directly, which does accept a per-subplot
        // List<Boolean>.
        //
        // That lower-level call builds a correct new chart internally,
        // but -- unlike GraphWidget.drawGraph() -- does NOT wire it into
        // GraphPanel's own visible chartPane itself; it leaves its
        // internal graphOn flag false and returns, expecting its caller
        // to flip that and swap in the new chart (see
        // GraphPanel.togglePlot()). GraphWidget.drawGraph() normally does
        // exactly that immediately after its own single-flag call to
        // drawGraphPanel(); since we call drawGraphPanel directly here
        // instead, going around GraphWidget.drawGraph(), both steps below
        // are required: setPlotSpec() first (so GraphWidget's own
        // isPlotEmpty() check -- which looks at GraphWidget's OWN
        // plotSpecs field, unrelated to GraphPanel's internal state --
        // sees a real plot and agrees to show GraphPanel at all), then
        // plotWidget.togglePlot() after the real draw call, to make
        // GraphPanel actually display the chart that was just built
        // instead of silently keeping whatever was on screen before.
        plotWidget.setPlotSpec(hazardSpec);

        plotWidget.getGraphPanel().drawGraphPanel(
                List.of(hazardSpec, fragilitySpec),
                List.of(true, true),    // X (IML) log on both subplots
                List.of(true, false),   // hazard log-Y, fragility linear-Y
                null, null);

        plotWidget.togglePlot();

        addVerticalMarkersAcrossSubplots(disaggBounds);
    }

    private DiscretizedFunc buildFunc(double[] x, double[] y, String name) {
        ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
        for (int i = 0; i < x.length; i++)
            func.set(x[i], y[i]);
        func.setName(name);
        return func;
    }

    /**
     * Adds a vertical {@link ValueMarker} at each of {@code imlValues}
     * to every subplot of {@link #plotWidget}'s combined plot --
     * {@code CombinedDomainXYPlot.addDomainMarker} only draws on the
     * subplot it's called on, so each subplot needs its own call with
     * the same X value to make the lines appear to span the whole
     * stacked figure.
     */
    private void addVerticalMarkersAcrossSubplots(double[] imlValues) {
        XYPlot plot = plotWidget.getGraphPanel().getPlot();

        List<XYPlot> subplots = (plot instanceof CombinedDomainXYPlot)
                ? ((CombinedDomainXYPlot) plot).getSubplots()
                : List.of(plot);

        for (double iml : imlValues) {
            ValueMarker marker = new ValueMarker(iml);
            marker.setPaint(Color.GRAY);
            marker.setStroke(new java.awt.BasicStroke(
                    1f, java.awt.BasicStroke.CAP_BUTT, java.awt.BasicStroke.JOIN_MITER,
                    1f, new float[] {6f, 4f}, 0f));
            for (XYPlot subplot : subplots)
                subplot.addDomainMarker(marker);
        }
    }

    /**
     * Ported from the original {@code PBRAppUI.PlotData}: linear
     * interpolation to find the IML at which {@code cdf} crosses each
     * of {@code markerValues}.
     */
    private static double[] imlsAtMarkers(double[] imls, double[] cdf, double[] markerValues) {
        double[] result = new double[markerValues.length];
        for (int m = 0; m < markerValues.length; m++)
            result[m] = imlAtProbability(imls, cdf, markerValues[m]);
        return result;
    }

    private static double imlAtProbability(double[] imls, double[] cdf, double targetProb) {
        int n = imls.length;
        int idx = java.util.Arrays.binarySearch(cdf, targetProb);
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

    private void logResult(String line) {
        resultsText.append(line + "\n");
    }

    public static void main(String[] args) {
        javax.swing.SwingUtilities.invokeLater(() -> new PBRAppUIv2().setVisible(true));
    }
}
