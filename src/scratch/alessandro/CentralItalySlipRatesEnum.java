/**
 * 
 */
package scratch.alessandro.logicTreeEnums;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.Ellsworth_B_WG02_MagAreaRel;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.HanksBakun2002_MagAreaRel;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.Shaw_2009_ModifiedMagAreaRel;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.GraphWindow;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.io.Files;

/**
 * * 
 * @author field
 *
 */
public enum CentralItalySlipRatesEnum {
	
//	UNIFORM_MEAN("UNIFORM_MEAN", "src/scratch/alessandro/data/uniform_deformation_model_0703.txt", 1.0) {},
//	
//	UNIFORM_LOWER("UNIFORM_LOWER", "src/scratch/alessandro/data/uniform_deformation_model_0703.txt", 0.76) {},
//
//	UNIFORM_UPPER("UNIFORM_UPPER", "src/scratch/alessandro/data/uniform_deformation_model_0703.txt", 1.28) {},
	
	// ok for the moment, but change upper and lower with the right values.
	NON_UNIFORM_MEAN("NON_UNIFORM_MEAN", "src/scratch/alessandro/data__CentralItaly/slip_rate_CentralItaly.txt", 1.0) {},
	
	NON_UNIFORM_LOWER("NON_UNIFORM_LOWER", "src/scratch/alessandro/data__CentralItaly/slip_rate_CentralItaly.txt", 0.76) {},

	NON_UNIFORM_UPPER("NON_UNIFORM_UPPER", "src/scratch/alessandro/data__CentralItaly/slip_rate_CentralItaly.txt", 1.28) {};


	private String name;
	double[] sectSlipRate, sectSlipRateStdDev;

	private CentralItalySlipRatesEnum(String name, String fileName, double scaleFactor) {
		this.name = name;

		try {
			File file = new File(fileName);
			List<String> fileLines = Files.readLines(file, Charset.defaultCharset());
			int numSections = fileLines.size();
			int sectIndex = 0;
			sectSlipRate = new double[numSections];
			sectSlipRateStdDev = new double[numSections];
			for (String line : fileLines) {
				//			System.out.println(line);
				line = line.trim();
				String[] split = line.split("\t");	// tab delimited
				Preconditions.checkState(split.length == 4, "Expected 4 items, got %s", split.length);
				//			System.out.println(split[0]+"\t"+split[1]+"\t"+split[2]);
				int testIndex = Integer.valueOf(split[0]);
				if(sectIndex != testIndex)
					throw new RuntimeException("Bad section index; "+sectIndex+" != "+testIndex);
				sectSlipRate[sectIndex] = scaleFactor*Double.valueOf(split[1]);
				double low95 = scaleFactor*Double.valueOf(split[2]);
				double upp95 = scaleFactor*Double.valueOf(split[3]);
				sectSlipRateStdDev[sectIndex] = (upp95-low95)/(2*1.96);	
				sectIndex+=1;
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	/**
	 * This returns the slip rate array (mm/yr)
	 * @return
	 */
	public  double[] getSlipRateArray() {
		return sectSlipRate;
	}

	/**
	 * This returns the slip rate standard deviation of the mean (mm/yr)
	 * @return
	 */
	public  double[] getSlipRateStdomArray() {
		return sectSlipRateStdDev;
	}


	public String getName() {
		return name;
	}

	public String getBranchLevelName() {
		return "Central Italy Slip Rates";
	}


	//public 
	public static void main(String[] args) throws IOException {


	}


}
