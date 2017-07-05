package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_SimAnalysisTools;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.ByteStreams;

public class RunawayDetect {
	
//	private static 

	public static void main(String[] args) throws ZipException, IOException {
		File zipFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2015_10_12-spontaneous-100yr-full_td-maxChar10.0/results.zip");
		
		ZipFile zip = new ZipFile(zipFile);

		List<List<ETAS_EqkRupture>> catalogs = Lists.newArrayList();
		
		ArrayList<? extends ZipEntry> entries = Collections.list(zip.entries());
		Collections.sort(entries, new Comparator<ZipEntry>() {

			@Override
			public int compare(ZipEntry o1, ZipEntry o2) {
				return o1.getName().compareTo(o2.getName());
			}
		});
		
		SimpleDateFormat df = new SimpleDateFormat("yyyy_MM_dd");
		
		System.setErr(new PrintStream(ByteStreams.nullOutputStream()));

		for (ZipEntry entry : entries) {
			if (!entry.isDirectory())
				continue;
			//			System.out.println(entry.getName());
			String subEntryName = entry.getName()+"simulatedEvents.txt";
			ZipEntry catEntry = zip.getEntry(subEntryName);
			if (catEntry == null)
				continue;
			//			System.out.println("Loading "+catEntry.getName());

			try {
				List<ETAS_EqkRupture> cat = ETAS_CatalogIO.loadCatalog(
						zip.getInputStream(catEntry), 0d, true);
				
				Preconditions.checkState(!cat.isEmpty());
				ETAS_EqkRupture last = cat.get(cat.size()-1);
				int year = last.getOriginTimeCal().get(Calendar.YEAR);
				List<ETAS_EqkRupture> supra = Lists.newArrayList();
				for (ETAS_EqkRupture rup : cat)
					if (rup.getFSSIndex() >= 0)
						supra.add(rup);
				int numSupra = supra.size();
				
				Map<String, Integer> dateCounts = Maps.newHashMap();
				for (ETAS_EqkRupture rup : supra) {
					String date = df.format(rup.getOriginTimeCal().getTime());
					Integer count = dateCounts.get(date);
					if (count == null)
						count = 0;
					count++;
					dateCounts.put(date, count);
				}
				
				String maxDate = null;
				int maxCount = 0;
				
				for (String date : dateCounts.keySet()) {
					int count = dateCounts.get(date);
					if (count > maxCount) {
						maxCount = count;
						maxDate = date;
					}
				}
				
				System.out.println(entry.getName()+"\tyearEnd="+year+"\tnumSupra="+numSupra
						+"\tmax # supra on a single day:\t"+maxCount+" on "+maxDate);

				catalogs.add(cat);
			} catch (Exception e) {
				//				ExceptionUtils.throwAsRuntimeException(e);
				System.out.println("Skipping catalog "+entry.getName()+": "+e.getMessage());
			}
			if (catalogs.size() % 1000 == 0)
				System.out.println("Loaded "+catalogs.size()+" catalogs (and counting)...");
		}

		zip.close();
	}

}
