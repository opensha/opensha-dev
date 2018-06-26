package scratch.kevin.simulators.synch;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JPanel;

import org.dom4j.DocumentException;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.iden.LogicalOrRupIden;
import org.opensha.sha.simulators.iden.RuptureIdentifier;

import com.google.common.base.Preconditions;

import scratch.kevin.simulators.PeriodicityPlotter;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimMarkovChainBuilder;
import scratch.kevin.util.MarkdownUtils;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class RSQSimOccupancyCopulaCalculator {

	public static void main(String[] args) throws IOException, DocumentException {
		File mainOutputDir = new File("/home/kevin/git/rsqsim-analysis/catalogs");
		File catalogsDir = new File("/data/kevin/simulators/catalogs");
		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(catalogsDir);
		
		List<String> parentNames = new ArrayList<>();
		
		parentNames.add("San Andreas (Carrizo) rev");
		parentNames.add("San Andreas (Cholame) rev");
		parentNames.add("San Andreas (Mojave S)");
		parentNames.add("San Andreas (Coachella) rev");
		parentNames.add("San Jacinto (Anza) rev");
		parentNames.add("Garlock (West)");
		
		double minMag = 7d;
		OccupancyCopulaCalculator.lowResTimeDelta = 5d; // years
		double skipYears = 5000;
		boolean middleSubSect = true; // else any
		
		int numCopulaBins = 50;
		
		List<RuptureIdentifier> faultIdens = RSQSimMarkovChainBuilder.getParentFaultIdens(catalog, middleSubSect, parentNames);
		List<RSQSimEvent> events = catalog.loader().minMag(minMag).skipYears(skipYears).matches(new LogicalOrRupIden(faultIdens)).load();
		
		File catOutDir = new File(mainOutputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catOutDir.exists() || catOutDir.mkdir());
		
		int dims = faultIdens.size();
		
		File catCopulaOutDir = new File(catOutDir, "occupancy_copula_m"+(float)minMag+"_"+dims+"D");
		Preconditions.checkState(catCopulaOutDir.exists() || catCopulaOutDir.mkdir());
		File resourcesDir = new File(catCopulaOutDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		List<String> lines = new ArrayList<>();
		
		// header
		lines.add("# "+dims+"-D M>="+(float)minMag+" Occupancy Copulas");
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		for (int i=0; i<faultIdens.size(); i++) {
			RuptureIdentifier iden1 = faultIdens.get(i);
			String name1 = iden1.getName();
			for (int j=i+1; j<faultIdens.size(); j++) {
				RuptureIdentifier iden2 = faultIdens.get(j);
				String name2 = iden2.getName();
				
				System.out.println("Doing "+name1+" and "+name2);
				
				lines.add("## "+name1+" vs "+name2);
				lines.add(topLink); lines.add("");
				lines.add("");
				
				OccupancyCopulaCalculator calc = new OccupancyCopulaCalculator(events, iden1, iden2);
//				calc.plotCopulaCombo();
				JPanel comboPlot = calc.getCopulaComboPlot(numCopulaBins);
				String prefix = "copula_"+PeriodicityPlotter.getFileSafeString(name1)
						+"_"+PeriodicityPlotter.getFileSafeString(name2);
				OccupancyCopulaCalculator.writePlotPNG(comboPlot, new File(resourcesDir, prefix+".png"));
				OccupancyCopulaCalculator.writePlotPDF(comboPlot, new File(resourcesDir, prefix+".pdf"));
				
				lines.add("![Occupancy Copula]("+resourcesDir.getName()+"/"+prefix+".png)");
				lines.add("");
			}
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");
		
		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, catCopulaOutDir);
		
		catalog.writeMarkdownSummary(catOutDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(mainOutputDir);
	}

}
