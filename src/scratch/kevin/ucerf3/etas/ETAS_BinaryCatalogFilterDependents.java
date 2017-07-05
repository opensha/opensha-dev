package scratch.kevin.ucerf3.etas;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.util.List;

import org.opensha.commons.util.ClassUtils;

import com.google.common.base.Preconditions;

import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO.BinarayCatalogsIterable;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_SimAnalysisTools;

public class ETAS_BinaryCatalogFilterDependents {

	public static void main(String[] args) {
		if (args.length < 2 || args.length > 3) {
			System.err.println("USAGE: "+ClassUtils.getClassNameWithoutPackage(ETAS_BinaryCatalogFilterDependents.class)
						+" <input-binary-file> <output-binary-file> [<trigger-id>]");
			System.exit(2);
		}
		
		File inputFile = new File(args[0]);
		Preconditions.checkArgument(inputFile.exists(), "Input file doesn't exist: %s", inputFile.getAbsolutePath());
		
		File outputFile = new File(args[1]);
//		Preconditions.checkArgument(outputFile.getParentFile().exists(),
//				"Output file cannot be created: %s", outputFile.getAbsolutePath());
		
		int triggerParentID = 0;
		if (args.length == 3)
			triggerParentID = Integer.parseInt(args[2]);
		
		try {
			BinarayCatalogsIterable iterable = ETAS_CatalogIO.getBinaryCatalogsIterable(inputFile, -1);
			int numCatalogs = iterable.getNumCatalogs();
			
			long inNumRups = 0;
			long outNumRups = 0;
			
			System.out.println("Processing "+numCatalogs+" catalogs, triggerParentID="+triggerParentID);
			
			DataOutputStream out = new DataOutputStream(new BufferedOutputStream(
					new FileOutputStream(outputFile), ETAS_CatalogIO.buffer_len));

			// write number of catalogs as int
			out.writeInt(numCatalogs); // will overwrite later
			
			int catalogIndex = 0;
			
			for (List<ETAS_EqkRupture> catalog : iterable) {
				inNumRups += catalog.size();
				List<ETAS_EqkRupture> children = ETAS_SimAnalysisTools.getChildrenFromCatalog(catalog, triggerParentID);
				outNumRups += children.size();
				
				if (catalogIndex++ % 1000 == 0)
					System.out.println("Processing catalog "+catalogIndex);
				
				ETAS_CatalogIO.writeCatalogBinary(out, children);
			}
			
			double percentChildren = 100d*(double)(outNumRups - inNumRups)/(double)inNumRups;
			
			System.out.println("Decendents represent "+(float)percentChildren+" % of input ruptures");
			
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.err.flush();
			System.exit(1);
		}
		System.exit(0);
	}

}
