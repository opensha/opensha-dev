package scratch.kevin.simulators.hazard;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.List;

import org.dom4j.Document;
import org.dom4j.Element;
import org.opensha.commons.param.ParamLinker;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.ParameterList;
import org.opensha.commons.param.event.ParameterChangeWarningListener;
import org.opensha.commons.param.impl.DoubleParameter;
import org.opensha.commons.param.impl.ParameterListParameter;
import org.opensha.commons.util.XMLUtils;
import org.opensha.sha.imr.AbstractIMR;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.MultiIMR_Averaged_AttenRel;
import org.opensha.sha.imr.mod.ModAttenRelRef;
import org.opensha.sha.imr.mod.ModAttenuationRelationship;
import org.opensha.sha.imr.mod.impl.FixedStdDevMod;

import com.google.common.base.Preconditions;

public class FixedStdDevNGAW2_GMPE extends MultiIMR_Averaged_AttenRel {
	
	private DoubleParameter stdDevParam;
	
	public FixedStdDevNGAW2_GMPE(ParameterChangeWarningListener listener) {
		// this constructor needed for loading from XML
		this(0.5);
	}

	public FixedStdDevNGAW2_GMPE(double stdDev) {
		super(getGMPEs());
		
		setParamDefaults();
		List<ModAttenuationRelationship> gmpes = (List<ModAttenuationRelationship>) this.getIMRs();
		
		stdDevParam = new DoubleParameter("Std. Dev.", 0.5);
		stdDevParam.setDefaultValue(0.5);
		stdDevParam.setValue(stdDev);
		
		for (ModAttenuationRelationship gmpe : gmpes) {
			ParameterListParameter pListParam =
					(ParameterListParameter) gmpe.getParameter("Modifier Params");
			Parameter<Double> oStdDev = pListParam.getParameter().getParameter(Double.class, FixedStdDevMod.STD_DEV_PARAM_NAME);
			new ParamLinker<Double>(stdDevParam, oStdDev);
			Preconditions.checkState(oStdDev.getValue() == stdDev);
		}
		
		otherParams.addParameter(stdDevParam);
	}
	
	private static List<ModAttenuationRelationship> getGMPEs() {
		List<ModAttenuationRelationship> gmpes = new ArrayList<>();
		
		gmpes.add(getGMPE(AttenRelRef.ASK_2014));
		gmpes.add(getGMPE(AttenRelRef.BSSA_2014));
		gmpes.add(getGMPE(AttenRelRef.CB_2014));
		gmpes.add(getGMPE(AttenRelRef.CY_2014));
		
		return gmpes;
	}
	
	private static ModAttenuationRelationship getGMPE(AttenRelRef ref) {
		ModAttenuationRelationship gmpe = new ModAttenuationRelationship(ref, ModAttenRelRef.FIXED_STD_DEV);
		gmpe.setParamDefaults();
		return gmpe;
	}
	
	public static void main(String[] args) throws InvocationTargetException {
		double stdDev = 0.23456;
		FixedStdDevNGAW2_GMPE origGMPE = new FixedStdDevNGAW2_GMPE(stdDev);
		
		System.out.println("**** Original GMPE ****");
		printStdDevDetails(origGMPE);
		
		Document doc = XMLUtils.createDocumentWithRoot();
		Element root = doc.getRootElement();
		origGMPE.toXMLMetadata(root);
		
		FixedStdDevNGAW2_GMPE newGMPE = (FixedStdDevNGAW2_GMPE) FixedStdDevNGAW2_GMPE.fromXMLMetadata(
				root.element(AbstractIMR.XML_METADATA_NAME), null);
		
		System.out.println("**** New GMPE ****");
		printStdDevDetails(newGMPE);
	}
	
	private static void printStdDevDetails(FixedStdDevNGAW2_GMPE gmpe) {
		System.out.println("Main param: "+gmpe.stdDevParam.getValue());
		System.out.println("Main getStdDev(): "+gmpe.getStdDev());
		for (ScalarIMR imr : gmpe.getIMRs()) {
			ModAttenuationRelationship modIMR = (ModAttenuationRelationship)imr;
			ParameterList modParamList = (ParameterList) modIMR.getParameter("Modifier Params").getValue();
			System.out.println(modIMR.getShortName());
			System.out.println("\tParam: "+modParamList.getParameter(FixedStdDevMod.STD_DEV_PARAM_NAME).getValue());
			System.out.println("\tgetStdDev(): "+imr.getStdDev());
		}
	}

}
