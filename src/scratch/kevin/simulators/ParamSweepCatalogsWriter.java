package scratch.kevin.simulators;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;

import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;

public class ParamSweepCatalogsWriter {

	public static void main(String[] args) throws IOException, DocumentException {
		File gitDir = new File("/home/scec-02/kmilner/git/rsqsim-analysis/2018_param_sweep");
		File catalogsDir = new File("/home/scec-00/gilchrij/RSQSim/CISM/paramSweep/paramSweep");
		File csvFile = new File("/home/scec-02/kmilner/simulators/catalogs/rsqsim_param_sweep_catalogs.csv");
		GregorianCalendar cal = RSQSimCatalog.cal(2018, 9, 27);
		String author = "Jaqcui Gilchrist";;
		FaultModels fm = FaultModels.FM3_1;
		DeformationModels dm = DeformationModels.GEOLOGIC;
		boolean replot = false;
		if (args.length == 1)
			replot = Boolean.parseBoolean(args[0]);
		double plotMinMag = 5d;
		
		System.out.println("Min Mag: "+plotMinMag);
		System.out.println("Replot? "+replot);
		
		Preconditions.checkState(gitDir.exists() || gitDir.mkdir());
		
		CSVFile<String> csv = CSVFile.readFile(csvFile, false);
		
		int dirCol = 0;
		int detailsCol = 2;
		int[] paramCols = IntStream.rangeClosed(3, 11).toArray();
		int baseRow = 1;
		Map<String, String> baseParams = null;
		List<String> paramNamesSorted = null;
		
		for (int row=baseRow; row<csv.getNumRows(); row++) {
			String dirName = csv.get(row, dirCol);
			String description = csv.get(row, detailsCol);
			if (row == baseRow) {
				baseParams = getParams(csv, row, paramCols);
				paramNamesSorted = new ArrayList<>(baseParams.keySet());
				Collections.sort(paramNamesSorted);
			} else {
				Preconditions.checkNotNull(baseParams);
				Map<String, String> params = getParams(csv, row, paramCols);
				List<String> diffParams = new ArrayList<>();
				for (String name : paramNamesSorted) {
					if (!params.get(name).equals(baseParams.get(name)))
						diffParams.add(name+"="+params.get(name));
				}
				Preconditions.checkState(!diffParams.isEmpty());
				description += ": "+Joiner.on(", ").join(diffParams);
			}
			System.out.println("Processing: "+dirName);
			System.out.println("\t"+description);
			
			File catalogDir = new File(catalogsDir, dirName);
			Preconditions.checkState(catalogDir.exists());
			
			Stopwatch watch = Stopwatch.createStarted();
			
			RSQSimCatalog catalog = new RSQSimCatalog(catalogDir, dirName, author, cal, description, fm, dm);
			try {
				double duration = catalog.getDurationYears();
				System.out.println("\tCatalog length: "+(float)duration+" years");
			} catch (IOException e) {
				System.out.println("IOException with catalog, skipping: "+e.getMessage());
				watch.stop();
				continue;
			}
			
			File catGitDir = new File(gitDir, catalog.getCatalogDir().getName());
			Preconditions.checkState(catGitDir.exists() || catGitDir.mkdir());
			catalog.writeMarkdownSummary(catGitDir, true, replot, plotMinMag);
			
			watch.stop();
			long secs = watch.elapsed(TimeUnit.SECONDS);
			double mins = secs/60d;
			System.out.println("\tTook "+secs+" seconds = "+(float)mins+" minutes");
		}
		
		RSQSimCatalog.writeCatalogsIndex(gitDir, true, csv.get(baseRow, dirCol));
	}
	
	private static Map<String, String> getParams(CSVFile<String> csv, int row, int[] paramCols) {
		Map<String, String> params = new HashMap<>();
		
		for (int col : paramCols) {
			String name = csv.get(0, col);
			String val = csv.get(row, col);
			params.put(name, val);
		}
		
		return params;
	}

}
