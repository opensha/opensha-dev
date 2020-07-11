package scratch.kevin.simulators.plots;

import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectionMapping;

import scratch.UCERF3.analysis.FaultBasedMapGen;

public class SectParticipationNucleationPlot extends AbstractPlot {
	
	private RSQSimSubSectionMapper mapper;
	private double[] minMags;
	
	private double overallMinMag;
	
	private List<? extends FaultSection> subSects;
	private double[][] particRates;
	private double[][] nuclRates;

	public SectParticipationNucleationPlot(RSQSimSubSectionMapper mapper, double... minMags) {
		this.mapper = mapper;
		this.minMags = minMags;
		
		overallMinMag = StatUtils.min(minMags);
		
		subSects = mapper.getSubSections();
		particRates = new double[minMags.length][subSects.size()];
		nuclRates = new double[minMags.length][subSects.size()];
	}

	@Override
	protected void doProcessEvent(SimulatorEvent e) {
		double mag = e.getMagnitude();
		if (mag < overallMinMag)
			return;
		List<List<SubSectionMapping>> mappings = mapper.getAllSubSectionMappings(e);
		for (List<SubSectionMapping> parentMappings : mappings) {
			for (SubSectionMapping mapping : parentMappings) {
				int index = mapping.getSubSect().getSectionId();
				for (int m=0; m<minMags.length; m++)
					if (mag >= minMags[m])
						particRates[m][index]++;
			}
		}
		double earliestTime = Double.POSITIVE_INFINITY;
		SimulatorElement hypoElem = null;
		for (EventRecord rec : e) {
			List<SimulatorElement> patches = rec.getElements();
			double[] patchTimes = rec.getElementTimeFirstSlips();
			for (int i=0; i<patches.size(); i++) {
				if (patchTimes[i] < earliestTime) {
					earliestTime = patchTimes[i];
					hypoElem = patches.get(i);
				}
			}
		}
		int hypoIndex = mapper.getMappedSection(hypoElem).getSectionId();
		for (int m=0; m<minMags.length; m++)
			if (mag >= minMags[m])
				nuclRates[m][hypoIndex]++;
	}

	@Override
	public void finalizePlot() throws IOException {
		double duration = getCurrentDurationYears();
		// normalize to rates
		for (int i=0; i<minMags.length; i++) {
			for (int j=0; j<subSects.size(); j++) {
				particRates[i][j] /= duration;
				nuclRates[i][j] /= duration;
			}
		}
		
		CPT ratioCPT = new CPT(0d, 1d, new Color(230, 230, 230),
				Color.RED.brighter(), Color.RED.darker().darker());
		ratioCPT.setNanColor(Color.GRAY);
		
		List<LocationList> faults = new ArrayList<>();
		for (FaultSection sect : subSects)
			faults.add(sect.getFaultTrace());
		
		for (int m=0; m<minMags.length; m++) {
			double[] particRates = this.particRates[m];
			double[] nuclRates = this.nuclRates[m];
			
			double minRate = 1d/duration;
			double maxRate = StatUtils.max(particRates);
			
			CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(
					Math.log10(minRate), Math.log10(maxRate));
			
			Region region = new CaliforniaRegions.RELM_TESTING();

			GMT_Map particMap = FaultBasedMapGen.buildMap(cpt, faults, FaultBasedMapGen.log10(particRates),
					null, 1d, region, false, "Log@-10@- M>="
							+optionalDigitDF.format(minMags[m])+" Participation Rate");
			GMT_Map nuclMap = FaultBasedMapGen.buildMap(cpt, faults, FaultBasedMapGen.log10(nuclRates),
					null, 1d, region, false, "Log@-10@- M>="
							+optionalDigitDF.format(minMags[m])+" Nucleation Rate");
			double[] ratios = new double[particRates.length];
			for (int i=0; i<particRates.length; i++)
				ratios[i] = nuclRates[i]/particRates[i];
//			GMT_Map ratioMap = FaultBasedMapGen.buildMap(ratioCPT, faults, FaultBasedMapGen.log10(ratios),
//					null, 1d, region, false, "Log@-10@- M>="
//							+optionalDigitDF.format(minMags[m])+" Nucleation/Participation Ratio");
			GMT_Map ratioMap = FaultBasedMapGen.buildMap(ratioCPT, faults, ratios,
					null, 1d, region, false,
					"M>="+optionalDigitDF.format(minMags[m])+" Nucleation/Participation Ratio");
			
			String prefix = getOutputPrefix()+"_m"+optionalDigitDF.format(minMags[m]);
			
//			System.out.println("Plotting");
			try {
				FaultBasedMapGen.plotMap(getOutputDir(), prefix+"_partic_map", false, particMap);
				FaultBasedMapGen.plotMap(getOutputDir(), prefix+"_nucl_map", false, nuclMap);
				FaultBasedMapGen.plotMap(getOutputDir(), prefix+"_nucl_partic_ratio_map", false, ratioMap);
			} catch (GMT_MapException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
		}
	}

	@Override
	public Collection<SimulatorElement> getApplicableElements() {
		return null;
	}

}
