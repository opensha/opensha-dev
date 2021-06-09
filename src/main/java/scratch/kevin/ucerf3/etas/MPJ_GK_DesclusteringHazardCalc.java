package scratch.kevin.ucerf3.etas;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FileUtils;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.AsyncPostBatchHook;
import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.ned.GK_Declustering.U3ETAS_SimulationAnalysis;

public class MPJ_GK_DesclusteringHazardCalc extends MPJTaskCalculator {
	
	private File outputDir;
	
	private GriddedRegion reg;
	
	private static final double[] periods = { 0, 0.2, 1, 5 };
	private static final AttenRelRef gmpeRef = AttenRelRef.CB_2014;
	private static final double duration = 50;
	
	private static final long random_seed = 123456789l;
	
	private ArrayList<ObsEqkRupList> catalogsList;
	private ArrayList<ObsEqkRupList> declusteredCatalogsList;
	private ArrayList<ObsEqkRupList> catalogsRandmizedList;
	
	private ExecutorService exec;
	
	private ZipOutputStream zip;
	
	public MPJ_GK_DesclusteringHazardCalc(CommandLine cmd) throws IOException, DocumentException {
		super(cmd);
		
		File catalogFile = new File(cmd.getOptionValue("catalog-file"));
		Preconditions.checkArgument(catalogFile.exists(),
				"Catalog file doesn't exist: %s", catalogFile.getAbsolutePath());
		
		File fssFile = new File(cmd.getOptionValue("fss-file"));
		Preconditions.checkArgument(fssFile.exists(),
				"FSS file doesn't exist: %s", fssFile.getAbsolutePath());
		
		double gridSpacing = Double.parseDouble(cmd.getOptionValue("grid-spacing"));
		
		outputDir = new File(cmd.getOptionValue("output-dir"));
		Preconditions.checkArgument(rank > 0 || outputDir.exists() || outputDir.mkdir(),
				"Output directory doesn't exist and can't be created: %s",
				outputDir.getAbsolutePath());
		
		reg = new CaliforniaRegions.RELM_TESTING_GRIDDED(gridSpacing);
		
		Random random = new Random(random_seed);
		
		debug("Loading catalogs...");
		catalogsList = U3ETAS_SimulationAnalysis.loadCatalogs(fssFile, catalogFile, 5d);
		debug("Declustering catalog...");
		declusteredCatalogsList = U3ETAS_SimulationAnalysis.getGK_DeclusteredCatalog(catalogsList);
		debug("Randomizing catalogs...");
		catalogsRandmizedList = U3ETAS_SimulationAnalysis.getRandomizedCatalogs(catalogsList, random);
		
		exec = Executors.newFixedThreadPool(getNumThreads());
		this.postBatchHook = new Hook();
		zip = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream(outputDir.getAbsolutePath()+".zip")));
		zip.setMethod(ZipOutputStream.DEFLATED);
	}

	@Override
	protected int getNumTasks() {
		return reg.getNodeCount();
	}

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		List<Future<?>> futures = new ArrayList<>();
		
		for (int index : batch)
			futures.add(exec.submit(new CalcTask(index)));
		
		for (Future<?> future : futures)
			future.get();
	}
	
	private File getOutputDirectory(int index) {
		Location loc = reg.getLocation(index);
		String name = "node_"+index+"_"+(float)loc.getLatitude()+"_"+(float)loc.getLongitude();
		return new File(MPJ_GK_DesclusteringHazardCalc.this.outputDir, name);
	}
	
	private class CalcTask implements Runnable {
		
		private int index;

		public CalcTask(int index) {
			this.index = index;
		}

		@Override
		public void run() {
			File outputDir = getOutputDirectory(index);
			Preconditions.checkState(outputDir.mkdir());
			
			ScalarIMR gmpe = gmpeRef.instance(null);
			
			Location loc = reg.getLocation(index);
			
			List<String> header = new ArrayList<>();
			header.add("x");
			header.add("Full Mean");
			header.add("Full Min");
			header.add("Full Max");
			header.add("Full Lower 95%");
			header.add("Full Upper 95%");
			header.add("Poisson");
			header.add("Declustered Mean");
			header.add("Declustered Min");
			header.add("Declustered Max");
			header.add("Declustered Lower 95%");
			header.add("Declustered Upper 95%");
			header.add("Randomized Mean");
			header.add("Randomized Min");
			header.add("Randomized Max");
			header.add("Randomized Lower 95%");
			header.add("Randomized Upper 95%");
			
			for (double period : periods) {
				String prefix;
				if (period == 0d) {
					gmpe.setIntensityMeasure(PGA_Param.NAME);
					prefix = "pga";
				} else {
					gmpe.setIntensityMeasure(SA_Param.NAME);
					SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), period);
					prefix = "sa_"+(float)period;
				}
				
				debug("Calculating index "+index+", period="+(float)period);
				
				UncertainArbDiscDataset[] datasetsArray1 = U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(
						catalogsList, loc, duration, period, false, gmpe);
				ArbitrarilyDiscretizedFunc poissonCurve = U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogsPoisson(
						catalogsList, loc, duration, period, gmpe);
				UncertainArbDiscDataset[] datasetsArray2 = U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(
						declusteredCatalogsList, loc, duration, period, false, gmpe);
				UncertainArbDiscDataset[] datasetsArray3 = U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(
						catalogsRandmizedList, loc, duration, period, false, gmpe);
				
				CSVFile<String> csv = new CSVFile<>(true);
				csv.addLine(header);
				
				for (int i=0; i<poissonCurve.size(); i++) {
					double x = poissonCurve.getX(i);
					List<String> line = new ArrayList<>();
					line.add(x+"");
					line.add(datasetsArray1[0].getY(i)+"");
					line.add(datasetsArray1[0].getLowerY(i)+"");
					line.add(datasetsArray1[0].getUpperY(i)+"");
					line.add(datasetsArray1[1].getLowerY(i)+"");
					line.add(datasetsArray1[1].getUpperY(i)+"");
					line.add(poissonCurve.getY(i)+"");
					line.add(datasetsArray2[0].getY(i)+"");
					line.add(datasetsArray2[0].getLowerY(i)+"");
					line.add(datasetsArray2[0].getUpperY(i)+"");
					line.add(datasetsArray2[1].getLowerY(i)+"");
					line.add(datasetsArray2[1].getUpperY(i)+"");
					line.add(datasetsArray3[0].getY(i)+"");
					line.add(datasetsArray3[0].getLowerY(i)+"");
					line.add(datasetsArray3[0].getUpperY(i)+"");
					line.add(datasetsArray3[1].getLowerY(i)+"");
					line.add(datasetsArray3[1].getUpperY(i)+"");
					csv.addLine(line);
				}
				
				try {
					csv.writeToFile(new File(outputDir, prefix+".csv"));
				} catch (IOException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
			}
			debug("done with "+index);
		}
		
	}
	
	private class Hook extends AsyncPostBatchHook {

		public Hook() {
			super(1);
		}

		@Override
		protected void batchProcessedAsync(int[] batch, int processIndex) {
			debug("zipping files for "+batch.length+" nodes from "+processIndex+".\t"+getCountsString());
			
			int BUFFER = 8192;
			byte data[] = new byte[BUFFER];
			for (int index : batch) {
				File outputDir = getOutputDirectory(index);
				
				try {
					zip.putNextEntry(new ZipEntry(outputDir.getName()+"/"));
					for (File file : outputDir.listFiles()) {
						if (!file.getName().endsWith(".csv"))
							continue;
						zip.putNextEntry(new ZipEntry(outputDir.getName()+"/"+file.getName()));
						FileInputStream fi = new FileInputStream(file);
						BufferedInputStream origin = new BufferedInputStream(fi, BUFFER);
						int count;
						while ( (count = origin.read(data, 0,
								BUFFER)) != -1) {
							zip.write(data, 0, count);
						}
						origin.close();
					}
					FileUtils.deleteRecursive(outputDir);
				} catch (IOException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
			}
		}
		
	}

	@Override
	protected void doFinalAssembly() throws Exception {
		exec.shutdown();
		((Hook)this.postBatchHook).shutdown();
		zip.close();
	}
	
	public static Options createOptions() {
		Options ops = MPJTaskCalculator.createOptions();
		
		Option catalogs = new Option("c", "catalog-file", true, "ETAS Catalogs Binary File");
		catalogs.setRequired(true);
		ops.addOption(catalogs);
		
		Option spacing = new Option("s", "grid-spacing", true, "Grid spacing in degrees");
		spacing.setRequired(true);
		ops.addOption(spacing);
		
		Option output = new Option("o", "output-dir", true, "Output directory");
		output.setRequired(true);
		ops.addOption(output);
		
		Option fss = new Option("f", "fss-file", true, "Fault System Solution file");
		fss.setRequired(true);
		ops.addOption(fss);
		
		return ops;
	}
	
	public static void main(String args[]) {
		args = MPJTaskCalculator.initMPJ(args);
		
		try {
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_GK_DesclusteringHazardCalc.class);
			
			MPJ_GK_DesclusteringHazardCalc driver = new MPJ_GK_DesclusteringHazardCalc(cmd);
			
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
