package scratch.kevin.ucerf3.eal.spatialCorr;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.StatUtils;
import org.dom4j.DocumentException;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.jfree.chart.ui.TextAnchor;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_WrapperFullParam;
import org.opensha.sha.imr.attenRelImpl.ngaw2.ScalarGroundMotion;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.OtherParams.StdDevTypeParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sra.calc.EALCalculator;
import org.opensha.sra.gui.portfolioeal.Asset;
import org.opensha.sra.gui.portfolioeal.Portfolio;
import org.opensha.sra.vulnerability.Vulnerability;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.utils.FaultSystemIO;

public class SpatiallyCorrelatedLossDebug {

	public static void main(String[] args) throws IOException, DocumentException,
	ClassNotFoundException, InstantiationException, IllegalAccessException {
		FaultSystemSolution fss = FaultSystemIO.loadSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		File fieldDir = new File("/home/kevin/OpenSHA/UCERF3/eal/random_fields/sa10_1km_800x800");
		File portfolioFile = new File("/home/kevin/OpenSHA/UCERF3/eal/"
				+ "Porter-09-Feb-2020-CEA-100-pct-procy-portfolio-wills2015.csv");
		
		// randomly draw this many ruptures
		int numRups = 10;
		// randomly draw this many assets
		int numAssets = 5;
		// will use this many fields (first N)
		int numFields = 100;
		// number of randomly sampled between-event values
		int numBetweens = 1000;
		
		File outputDir = new File("/home/kevin/git/misc-research/loss_rand_fields");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		outputDir = new File(outputDir, numBetweens+"_betweens");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		FaultSystemRupSet rupSet = fss.getRupSet();
		Random r = new Random(rupSet.getNumRuptures());
		
		NormalDistribution normDist = new NormalDistribution(new Well19937c(rupSet.getNumRuptures()), 0d, 1d);
		double[] betweens = normDist.sample(numBetweens);
		
		System.out.println("mean between-event: "+(float)StatUtils.mean(betweens));
		
		List<String> lines = new ArrayList<>();
		lines.add("# Random Field Loss Debug");
		lines.add("");
		lines.add("Random field dubugging with the following parameters:");
		lines.add("");
		lines.add("* "+numRups+" randomly chosen ruptures");
		lines.add("* "+numAssets+" randomly chosen assets");
		lines.add("* "+numFields+" random fields");
		lines.add("* "+numBetweens+" randomly sampled (from a normal distribution) between-event values");
		lines.add("");
		lines.add("The mean of all "+numBetweens+" between-event values is "+(float)StatUtils.mean(betweens));
		lines.add("");
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		RandomFieldLoader[] fields = new RandomFieldLoader[numFields];
		DecimalFormat fieldDF = new DecimalFormat("000");
		System.out.println("Loading "+numFields+" random fields...");
		for (int f=0; f<numFields; f++)
			fields[f] = RandomFieldLoader.load(new File(fieldDir,
					"800x800SA10_"+fieldDF.format(f+1)+".csv"), 1d);
		System.out.println("Done loading");
		
		ScalarIMR gmpe = AttenRelRef.ASK_2014.instance(null);
		gmpe.setParamDefaults();
		
		Portfolio portfolio = Portfolio.createPortfolio(portfolioFile);
		List<Asset> assets = portfolio.getAssetList();
		Collections.shuffle(assets, r);
		assets = assets.subList(0, numAssets);
		
		SpatiallyCorrelatedLossCalc.D = false;
		RandomFieldLoader.D = false;
		
		for (int i=0; i<numRups; i++) {
			int rupIndex = r.nextInt(rupSet.getNumRuptures());
			
			double mag = rupSet.getMagForRup(rupIndex);
			double rate = rupSet.getAveRakeForRup(rupIndex);
			
			lines.add("## Rupture "+rupIndex+", M"+(float)mag);
			lines.add(topLink); lines.add("");
			
			System.out.println("Calculating for rupture "+rupIndex+" with M="+(float)mag);
			RuptureSurface surf = rupSet.getSurfaceForRupture(rupIndex, 1d);
			
			EqkRupture rup = new EqkRupture(mag, rate, surf, null);
			
			Location centroid = SpatiallyCorrelatedLossCalc.calcRupCentroid(surf);
			
			for (int a=0; a<assets.size(); a++) {
				Asset asset = assets.get(a);
				
				if (numAssets > 1) {
					lines.add("### Rupture "+rupIndex+", Asset"+a);
					lines.add(topLink); lines.add("");
				}
				
//				gmpe.setParamDefaults();
				Site site = new Site();
				for (Parameter<?> param : gmpe.getSiteParams())
					site.addParameter((Parameter<?>)param.clone());
				asset.siteSetup(site);
				site = asset.getSite();
				
				Vulnerability vulnModel = asset.getVulnModel();
				
				Preconditions.checkNotNull(vulnModel, "Vulnerability model is null!");
				String imt = vulnModel.getIMT();
				double imls[] = vulnModel.getIMLValues();
				
				gmpe.setIntensityMeasure(imt);
				if (imt.equals(SA_Param.NAME))
					SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), vulnModel.getPeriod());
				System.out.println("\t\tAsset with value="+(float)asset.getValue());
				gmpe.setAll(rup, site, gmpe.getIntensityMeasure());
				double expectedVs30 = (double)asset.getParameterList().getParameter("Vs30").getValue();
				double siteVs30 = (Double)site.getParameter(Vs30_Param.NAME).getValue();
				double gmpeVs30 = (Double)site.getParameter(Vs30_Param.NAME).getValue();
				Preconditions.checkState(expectedVs30 == siteVs30, "vs30 should be %s, site has %s",
						expectedVs30, siteVs30);
				Preconditions.checkState(expectedVs30 == gmpeVs30, "vs30 should be %s, GMPE has %s",
						expectedVs30, gmpeVs30);
				double meanIML, stdDev, phi, tau;
				if (gmpe instanceof NGAW2_WrapperFullParam) {
					ScalarGroundMotion gm = ((NGAW2_WrapperFullParam)gmpe).getGroundMotion();
					meanIML = gm.mean();
					phi = gm.phi();
					tau = gm.tau();
					stdDev = gm.stdDev();
				} else {
					meanIML = gmpe.getMean();
					stdDev = gmpe.getStdDev();
					StdDevTypeParam type = (StdDevTypeParam)(gmpe.getParameter(StdDevTypeParam.NAME));
					type.setValue(StdDevTypeParam.STD_DEV_TYPE_INTER);
					tau = gmpe.getStdDev();
					type.setValue(StdDevTypeParam.STD_DEV_TYPE_INTRA);
					phi = gmpe.getStdDev();
				}
				
				lines.add("Asset parameters:");
				lines.add("");
				lines.add("* value: "+asset.getValue());
				lines.add("* J-B distance to rupture: "+rup.getRuptureSurface().getDistanceJB(site.getLocation()));
				lines.add("* Vulnerability function: "+asset.getVulnModelName());
				lines.add("* GMPE ln mean: "+(float)meanIML);
				lines.add("* GMPE phi/tau/total std devs: "+(float)phi+"\t"+(float)tau+"\t"+(float)stdDev);
				lines.add("");
				
				TableBuilder lossTable = MarkdownUtils.tableBuilder();
				lossTable.addLine("Description", "Loss Value");
				
				ArbitrarilyDiscretizedFunc logHazFunction = new ArbitrarilyDiscretizedFunc();
				ArbitrarilyDiscretizedFunc hazFunction = new ArbitrarilyDiscretizedFunc();

				// Setup for the HazardCurveCalculator
				HazardCurveCalculator calc = new HazardCurveCalculator();

				// Take the log of the x values of the hazard function
				// Used to make calculations
				for( int l = 0; l < imls.length; ++l ) {
//					System.out.println(i+". "+imls[i]+", log: "+Math.log(imls[i]));
					logHazFunction.set( Math.log( imls[l] ), 1 );
					hazFunction.set(imls[l], 1);
				}
				
				EALCalculator currentCalc = new EALCalculator(hazFunction,
						vulnModel.getVulnerabilityFunc(), asset.getValue());
				
				// calc deterministic hazard curve
				calc.getHazardCurve(logHazFunction, site, gmpe, rup);
				
				Preconditions.checkState(imls.length == logHazFunction.size());

				// populate the linear func with the y values
				for (int l=0; l<logHazFunction.size(); l++)
					hazFunction.set(l, logHazFunction.getY(l));

				currentCalc.setMAFE(hazFunction);
				double rupEL = currentCalc.computeEAL();
				lossTable.addLine("Mean Loss (GMPE)", (float)rupEL);
				
				NormalDistribution gmpeDist = new NormalDistribution(meanIML, stdDev);

				EvenlyDiscretizedFunc gmpeHist = new EvenlyDiscretizedFunc(meanIML-4*stdDev, meanIML+4*stdDev, 100);
				EvenlyDiscretizedFunc fieldHist = gmpeHist.deepClone();
				for (int j=0; j<gmpeHist.size(); j++) {
					double x = gmpeHist.getX(j);
					double x0 = x - 0.5*gmpeHist.getDelta();
					double x1 = x + 0.5*gmpeHist.getDelta();
					gmpeHist.set(j, gmpeDist.probability(x0, x1));
				}
				
				System.out.println("\t\t\tmean loss="+(float)rupEL);
				
				System.out.println("\t\t\tLn x-values:\n\t\t\t\t"+Joiner.on(",").join(logHazFunction.xValues()));
				System.out.println("\t\t\tProb of exceed:\n\t\t\t\t"+Joiner.on(",").join(hazFunction.yValues()));
				
				List<Asset> myAssets = new ArrayList<>();
				myAssets.add(asset);
				Table<Double, RandomFieldLoader, Double> vals = SpatiallyCorrelatedLossCalc.calcSpatiallyCorrelatedLoss(
									gmpe, myAssets, rup, centroid, betweens, fields);
				double[] randLossVals = Doubles.toArray(vals.values());
				double fieldsMean = StatUtils.mean(randLossVals);
				System.out.println("\tmean of "+randLossVals.length+" randoms: "+(float)fieldsMean);
				lossTable.addLine("Mean Loss ("+randLossVals.length+" random fields)", (float)fieldsMean);
				double[] logVals = new double[randLossVals.length];
				for (int j=0; j<logVals.length; j++)
					logVals[j] = Math.log(randLossVals[j]);
				double logMean = StatUtils.mean(logVals);
				System.out.println("\tln mean of randoms: "+(float)logMean);
				
				// calculate exceedance curve from random fields
				DiscretizedFunc randFieldExceed = logHazFunction.deepClone();
				for (int j=0; j<randFieldExceed.size(); j++)
					randFieldExceed.set(j, 0d);
				// also figure out what the average random field value is here
				double[] randFieldVals = new double[fields.length];
				List<Double> randFieldIMLs = new ArrayList<>();
				
				ArbitrarilyDiscretizedFunc fieldCumulativeLosses = new ArbitrarilyDiscretizedFunc();
				for (int j=0; j<hazFunction.size(); j++)
					fieldCumulativeLosses.set(hazFunction.getX(j), 0d);
				
				for (int f=0; f<fields.length; f++) {
					randFieldVals[f] = fields[f].getValue(site.getLocation(), centroid);
					for (double between : betweens) {
						// ground motion considering the between-event term
						double tauGM = meanIML + between*tau;
						// ground motion considering the randomly sampled within-event term (and also tau)
						double phiGM = tauGM + randFieldVals[f]*phi;
						randFieldIMLs.add(phiGM);
						fieldHist.add(fieldHist.getClosestXIndex(phiGM), 1d);
						Double randLoss = vals.get(between, fields[f]);
						for (int j=0; j<randFieldExceed.size(); j++) {
							if (phiGM >= randFieldExceed.getX(j))
								randFieldExceed.set(j, randFieldExceed.getY(j)+1d);
							if (phiGM <= randFieldExceed.getX(j))
								fieldCumulativeLosses.set(j, fieldCumulativeLosses.getY(j)+randLoss);
						}
					}
				}
				fieldHist.scale(1d/randFieldIMLs.size());
				randFieldExceed.scale(1d/randFieldIMLs.size());
				fieldCumulativeLosses.scale(1d/randFieldIMLs.size());
				System.out.println("\t\t\tmean random field value: "+(float)StatUtils.mean(randFieldVals));
				System.out.println("\t\t\tmean IML from random field ground motions: "
						+(float)StatUtils.mean(Doubles.toArray(randFieldIMLs)));
				System.out.println("\t\t\texceedance curve calculated from random fields:\n\t\t\t\t"
						+Joiner.on(",").join(randFieldExceed.yValues()));
				
				// now calculate using the random field-based exceedence curve

				// populate the linear func with the y values
				DiscretizedFunc randFieldHazFunc = hazFunction.deepClone();
				for (int l=0; l<randFieldHazFunc.size(); l++)
					randFieldHazFunc.set(l, randFieldExceed.getY(l));

				currentCalc.setMAFE(randFieldHazFunc);
				double randEL = currentCalc.computeEAL();
				lossTable.addLine("Mean Loss using rand-field exceed probs", (float)randEL);
				lines.add("**Loss Values**");
				lines.add("");
				lines.addAll(lossTable.build());
				lines.add("");
				System.out.println("\t\t\tmean using the rand field exceed curve: "+(float)randEL);
				
//				TableBuilder xValsTable = MarkdownUtils.tableBuilder();
//				xValsTable.addLine("IML", "Ln IML", "Mean Damage Factor", "Exceed Prob (from GMPE)",
//						"Loss for GMPE exceed prob", "Exceed Prob (from random fields)", "Loss for RF exceed prob");
//				DiscretizedFunc vulnFunc = vulnModel.getVulnerabilityFunc();
//				for (int l=0; l<hazFunction.size(); l++)
//					xValsTable.addLine((float)hazFunction.getX(l), (float)logHazFunction.getX(l),
//							(float)vulnModel.getMeanDamageFactor(hazFunction.getX(l)),
//							(float)hazFunction.getY(l), (float)getLossFromBin(asset, hazFunction, l),
//							(float)randFieldExceed.getY(l), (float)getLossFromBin(asset, randFieldExceed, l));
//				lines.add("**Exceedance Probs**");
//				lines.add("");
//				lines.addAll(xValsTable.build());
//				lines.add("");
				
				TableBuilder plotsTable = MarkdownUtils.tableBuilder();
				plotsTable.initNewLine();
				
				String commonPrefix = "rup"+i+"_asset"+a;
				
				HeadlessGraphPanel gp = new HeadlessGraphPanel();
				gp.setTickLabelFontSize(18);
				gp.setAxisLabelFontSize(24);
				gp.setPlotLabelFontSize(24);
				gp.setLegendFontSize(22);
				gp.setBackgroundColor(Color.WHITE);
				
				for (boolean logY : new boolean[] {false, true}) {
					List<XY_DataSet> funcs = new ArrayList<>();
					List<PlotCurveCharacterstics> chars = new ArrayList<>();
					
					fieldHist.setName("random field GM dist");
					funcs.add(fieldHist);
					chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLUE));
					
					gmpeHist.setName("GMPE log-normal dist");
					funcs.add(gmpeHist);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GREEN));
					
					PlotSpec spec = new PlotSpec(funcs, chars, "GMPE Distributions", "Ln IML", "Density");
					spec.setLegendVisible(true);
					
					Range xRange = new Range(gmpeHist.getMinX(), gmpeHist.getMaxX());
					Range yRange;
					String logAdd = "";
					double maxY = Math.max(1.2*gmpeHist.getMaxY(), 1.05*fieldHist.getMaxY());
					if (logY) {
						double minY = gmpeHist.getMinY();
						if (minY == 0)
							minY = 1e-6;
						else
							minY = Math.pow(10, Math.floor(Math.log10(minY)));
						yRange = new Range(minY, maxY);
						logAdd = "_log";
					} else {
						yRange = new Range(0d, maxY);
					}
					
					addMeanSDAnns(Doubles.toArray(randFieldIMLs), spec, xRange, yRange, logY);
					
					gp.drawGraphPanel(spec, false, logY, xRange, yRange);
					
					File file = new File(resourcesDir, commonPrefix+"_gmpe_dists"+logAdd);
					gp.getChartPanel().setSize(800, 600);
					gp.saveAsPNG(file.getAbsolutePath()+".png");
					gp.saveAsPDF(file.getAbsolutePath()+".pdf");
					
					plotsTable.addColumn("![plot](resources/"+file.getName()+".png)");
					
					funcs = new ArrayList<>();
					chars = new ArrayList<>();
					
					double maxLoss = 2d*Math.max(rupEL, randEL);
					maxLoss = Math.max(maxLoss, 1.1*StatUtils.max(randLossVals));
					System.out.println("Histogram max loss: "+maxLoss);
					EvenlyDiscretizedFunc lossHist = HistogramFunction.getEncompassingHistogram(
							0d, maxLoss, maxLoss/100d);
					for (double lossVal : randLossVals)
						lossHist.add(lossHist.getClosestXIndex(lossVal), 1d);
					
					lossHist.setName("Random Field Losses");
					funcs.add(lossHist);
					chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));
					
					xRange = new Range(0d, maxLoss);
					
					yRange = new Range(0d, lossHist.getMaxY()*1.05d);
					if (logY)
						yRange = new Range(0.9, yRange.getUpperBound());
					
					DefaultXY_DataSet meanLossXY = new DefaultXY_DataSet();
					meanLossXY.set(rupEL, yRange.getLowerBound());
					meanLossXY.set(rupEL, yRange.getUpperBound());
					
					meanLossXY.setName("Mean Loss");
					funcs.add(meanLossXY);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GREEN));
					
					DefaultXY_DataSet fieldLossXY = new DefaultXY_DataSet();
					fieldLossXY.set(fieldsMean, yRange.getLowerBound());
					fieldLossXY.set(fieldsMean, yRange.getUpperBound());
					
					fieldLossXY.setName("Mean Field Loss");
					funcs.add(fieldLossXY);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.RED));
					
					spec = new PlotSpec(funcs, chars, "Random Field Loss Samples", "Loss", "Count");
					spec.setLegendVisible(true);
					
					gp.drawGraphPanel(spec, false, logY, xRange, yRange);
					
					file = new File(resourcesDir, commonPrefix+"_field_losses"+logAdd);
					gp.getChartPanel().setSize(800, 600);
					gp.saveAsPNG(file.getAbsolutePath()+".png");
					gp.saveAsPDF(file.getAbsolutePath()+".pdf");
					
					plotsTable.addColumn("![plot](resources/"+file.getName()+".png)");
				}
				
				// between-event histogram
				
				List<XY_DataSet> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				EvenlyDiscretizedFunc betweenHist = new EvenlyDiscretizedFunc(-4d, 4d, 100);
				
				for (double between : betweens)
					betweenHist.add(betweenHist.getClosestXIndex(between), 1d);
				
				funcs.add(betweenHist);
				chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 3f, Color.MAGENTA));
				
				funcs.add(stdNorm(betweenHist));
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GRAY));
				
				PlotSpec spec = new PlotSpec(funcs, chars, "Between-Event Samples", "Residual", "Count");
				
				Range xRange = new Range(betweenHist.getMinX(), betweenHist.getMaxX());
				Range yRange = new Range(0d, 1.05d*betweenHist.getMaxY());
				
				addMeanSDAnns(betweens, spec, xRange, yRange, false);
				
				gp.drawGraphPanel(spec, false, false, xRange, yRange);
				
				File file = new File(resourcesDir, commonPrefix+"_betweens");
				gp.getChartPanel().setSize(800, 600);
				gp.saveAsPNG(file.getAbsolutePath()+".png");
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
				
				plotsTable.addColumn("![plot](resources/"+file.getName()+".png)");
				
				// within-event histogram
				
				funcs = new ArrayList<>();
				chars = new ArrayList<>();
				
				EvenlyDiscretizedFunc withinTotalHist = new EvenlyDiscretizedFunc(-4d, 4d, 100);
				
				for (double within : randFieldVals)
					withinTotalHist.add(withinTotalHist.getClosestXIndex(within), 1d);
				
				funcs.add(withinTotalHist);
				chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 3f, Color.RED));
				
				funcs.add(stdNorm(withinTotalHist));
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GRAY));
				
				spec = new PlotSpec(funcs, chars, "Within-Event Samples", "Residual", "Count");
				
				xRange = new Range(withinTotalHist.getMinX(), withinTotalHist.getMaxX());
				yRange = new Range(0d, 1.05d*withinTotalHist.getMaxY());
				
				addMeanSDAnns(randFieldVals, spec, xRange, yRange, false);
				
				gp.drawGraphPanel(spec, false, false, xRange, yRange);
				
				file = new File(resourcesDir, commonPrefix+"_withins");
				gp.getChartPanel().setSize(800, 600);
				gp.saveAsPNG(file.getAbsolutePath()+".png");
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
				
				plotsTable.addColumn("![plot](resources/"+file.getName()+".png)");

				funcs = new ArrayList<>();
				chars = new ArrayList<>();
				
				randFieldHazFunc.setName("RF exceed probs");
				funcs.add(randFieldHazFunc);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE));
				
				hazFunction.setName("GMPE exceed probs");
				funcs.add(hazFunction);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GREEN));
				
				DiscretizedFunc vulnFunc = vulnModel.getVulnerabilityFunc();
				vulnFunc.setName("Vuln Func");
				funcs.add(vulnFunc);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.GRAY));
				
				spec = new PlotSpec(funcs, chars, "GMPE Exceedance Probs", "IML", "Exceedance Prob");
				spec.setLegendVisible(true);
				
				double minX = Math.pow(10, Math.floor(Math.log10(hazFunction.getFirstInterpolatedX(0.999))));
				double maxX = Math.pow(10, Math.ceil(Math.log10(hazFunction.getFirstInterpolatedX(0.001))));
				minX = Math.max(minX, 1e-5);
				maxX = Math.min(maxX, 1e1);
				xRange = new Range(minX, maxX);
				yRange = new Range(0d, 1.02d);
				
				gp.drawGraphPanel(spec, true, false, xRange, yRange);
				
				file = new File(resourcesDir, commonPrefix+"_gmpe_exceeds");
				gp.getChartPanel().setSize(800, 600);
				gp.saveAsPNG(file.getAbsolutePath()+".png");
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
				
				plotsTable.addColumn("![plot](resources/"+file.getName()+".png)");
				
				// cumulative losses
				
				funcs = new ArrayList<>();
				chars = new ArrayList<>();
				
				double value = asset.getValue();
				
				DiscretizedFunc cmlLossesRFExceed = getCumulativeLosses(value, vulnFunc, randFieldHazFunc);
				
				cmlLossesRFExceed.setName("RF exceed probs");
				funcs.add(cmlLossesRFExceed);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE));
				
				DiscretizedFunc cmlLossesGMPEExceed = getCumulativeLosses(value, vulnFunc, hazFunction);
				
				cmlLossesGMPEExceed.setName("GMPE exceed probs");
				funcs.add(cmlLossesGMPEExceed);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GREEN));
				
				fieldCumulativeLosses.setName("Raw RF losses");
				funcs.add(fieldCumulativeLosses);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.RED));
				
				spec = new PlotSpec(funcs, chars, "Cumulative Losses", "IML", "Cumulative Loss");
				spec.setLegendVisible(true);
				
				yRange = new Range(0d, 1.1*Math.max(randEL, Math.max(rupEL, fieldsMean)));
				
				gp.drawGraphPanel(spec, true, false, xRange, yRange);
				
				file = new File(resourcesDir, commonPrefix+"_cumulative_losses");
				gp.getChartPanel().setSize(800, 600);
				gp.saveAsPNG(file.getAbsolutePath()+".png");
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
				
				plotsTable.addColumn("![plot](resources/"+file.getName()+".png)");
				
				plotsTable.finalizeLine();
				
				lines.addAll(plotsTable.wrap(2, 0).build());
				lines.add("");
				siteVs30 = (Double)site.getParameter(Vs30_Param.NAME).getValue();
				gmpeVs30 = (Double)site.getParameter(Vs30_Param.NAME).getValue();
				Preconditions.checkState(expectedVs30 == siteVs30, "vs30 should be %s, site has %s",
						expectedVs30, siteVs30);
				Preconditions.checkState(expectedVs30 == gmpeVs30, "vs30 should be %s, GMPE has %s",
						expectedVs30, gmpeVs30);
			}
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2, 4));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		System.out.print("Writing markdown...");
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
		System.out.println("DONE");
	}
	
	private static double getLossFromBin(Asset asset, DiscretizedFunc hazFunc, int index) {
		Vulnerability vuln;
		try {
			vuln = asset.getVulnModel();
		} catch (ClassNotFoundException | InstantiationException | IllegalAccessException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		return asset.getValue()*vuln.getMeanDamageAtExceedProb(hazFunc.getX(index), hazFunc.getY(index));
	}
	
	private static boolean super_sample_cumulative = false;
	
	private static DiscretizedFunc getCumulativeLosses(double value, DiscretizedFunc vulnFunc,
			DiscretizedFunc hazFunc) {
		if (super_sample_cumulative) {
			// resample
			double logMinX = Math.log(vulnFunc.getMinX());
			double logMaxX = Math.log(vulnFunc.getMaxX());
			EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(logMinX, logMaxX, 1000);
			
			DiscretizedFunc sampledDF = new ArbitrarilyDiscretizedFunc();
			DiscretizedFunc sampledHF = new ArbitrarilyDiscretizedFunc();
			for (int i=0; i<func.size(); i++) {
				double lnX = func.getX(i);
				double x = Math.exp(lnX);
				if (i == 0) {
					sampledDF.set(x, vulnFunc.getY(0));
					sampledHF.set(x, hazFunc.getY(0));
				} else if (i == func.size()-1) {
					sampledDF.set(x, vulnFunc.getY(vulnFunc.size()-1));
					sampledHF.set(x, hazFunc.getY(hazFunc.size()-1));
				} else {
					sampledDF.set(x, vulnFunc.getInterpolatedY(x));
					sampledHF.set(x, hazFunc.getInterpolatedY_inLogXDomain(x));
				}
			}
			vulnFunc = sampledDF;
			hazFunc = sampledHF;
		}
		
		DiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
		DiscretizedFunc modHazFunc = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<hazFunc.size(); i++)
			modHazFunc.set(hazFunc.getX(i), 0d);
		
		EALCalculator calc = new EALCalculator(hazFunc, vulnFunc, value);
		
		for (int i=0; i<hazFunc.size(); i++) {
			modHazFunc.set(i, hazFunc.getY(i));
			calc.setMAFE(modHazFunc);
			ret.set(hazFunc.getX(i), calc.computeEAL());
		}
		
		return ret;
	}
	
	private static double getLossFromBin(EALCalculator calc, DiscretizedFunc hazFunc, int index) {
		// first calculate it without that bin
		DiscretizedFunc newHazFunc = hazFunc.deepClone();
		for (int i=0; i<hazFunc.size(); i++) {
			if (i < index)
				newHazFunc.set(i, hazFunc.getY(i));
			else
				newHazFunc.set(i, 0d);
		}
			
		double lossBefore;
		if (index == 0) {
			lossBefore = 0d;
		} else {
			calc.setMAFE(newHazFunc);
			lossBefore = calc.computeEAL();
		}
		
		// now add that bin and calculate loss
		newHazFunc.set(index, hazFunc.getY(index));
		calc.setMAFE(newHazFunc);
		
		return calc.computeEAL() - lossBefore;
	}
	
	private static EvenlyDiscretizedFunc stdNorm(EvenlyDiscretizedFunc refHist) {
		EvenlyDiscretizedFunc ret = refHist.deepClone();
		double scale = refHist.calcSumOfY_Vals();
		
		NormalDistribution normDist = new NormalDistribution();
		
		double halfDelta = 0.5*ret.getDelta();
		for (int i=0; i<ret.size(); i++) {
			double x = ret.getX(i);
			double x0 = x-halfDelta;
			double x1 = x+halfDelta;
			ret.set(i, scale*normDist.probability(x0, x1));
		}
		return ret;
	}
	
	private static void addMeanSDAnns(double[] values, PlotSpec spec, Range xRange, Range yRange,
			boolean logY) {
		double mean = StatUtils.mean(values);
		double sd = Math.sqrt(StatUtils.variance(values));
		
		double x = 0.975*(xRange.getUpperBound() - xRange.getLowerBound()) + xRange.getLowerBound();
		double y1, y2, y3;
		if (logY) {
			double logLower = Math.log10(yRange.getLowerBound());
			double logUpper = Math.log10(yRange.getUpperBound());
			y1 = 0.975*(logUpper - logLower) + logLower;
			y1 = Math.pow(10, y1);
			y2 = 0.92*(logUpper - logLower) + logLower;
			y2 = Math.pow(10, y2);
			y3 = 0.865*(logUpper - logLower) + logLower;
			y3 = Math.pow(10, y3);
		} else {
			y1 = 0.975*(yRange.getUpperBound() - yRange.getLowerBound()) + yRange.getLowerBound();
			y2 = 0.92*(yRange.getUpperBound() - yRange.getLowerBound()) + yRange.getLowerBound();
			y3 = 0.865*(yRange.getUpperBound() - yRange.getLowerBound()) + yRange.getLowerBound();
		}
		
		DecimalFormat df = new DecimalFormat("0.00");
		
		Font font = new Font(Font.SANS_SERIF, Font.BOLD, 22);
		XYTextAnnotation meanAnn = new XYTextAnnotation("mean="+df.format(mean), x, y1);
		meanAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
		meanAnn.setFont(font);
		spec.addPlotAnnotation(meanAnn);
		XYTextAnnotation sdAnn = new XYTextAnnotation("sd="+df.format(sd), x, y2);
		sdAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
		sdAnn.setFont(font);
		spec.addPlotAnnotation(sdAnn);
		XYTextAnnotation nAnn = new XYTextAnnotation("N="+values.length, x, y3);
		nAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
		nAnn.setFont(font);
		spec.addPlotAnnotation(nAnn);
	}
	
	private static double calcLossAlt(double value, double sampledGM, DiscretizedFunc vulnFunc) {
		double df = vulnFunc.getInterpolatedY(sampledGM);
		
		return df*value;
	}

}
