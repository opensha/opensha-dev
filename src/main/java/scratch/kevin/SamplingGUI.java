package scratch.kevin;

import java.awt.BorderLayout;
import java.awt.Color;
import java.util.List;

import javax.swing.JFrame;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.gui.plot.GraphWidget;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotElement;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.param.ParameterList;
import org.opensha.commons.param.editor.impl.GriddedParameterListEditor;
import org.opensha.commons.param.event.ParameterChangeEvent;
import org.opensha.commons.param.event.ParameterChangeListener;
import org.opensha.commons.param.impl.ButtonParameter;
import org.opensha.commons.param.impl.DoubleParameter;

import com.google.common.collect.Lists;

public class SamplingGUI extends JFrame implements ParameterChangeListener {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private DoubleParameter toneFreqParam;
	private DoubleParameter samplingFreqParam;
	private DoubleParameter fractOffsetParam;
	private ButtonParameter randOffsetParam;
	
	private static final int wavelengths_to_plot = 10;
	
	private GraphWidget gw;
	
	public SamplingGUI() {
		setLayout(new BorderLayout());
		this.setSize(1000, 800);
		this.setDefaultCloseOperation(EXIT_ON_CLOSE);
		
		ParameterList params = new ParameterList();
		
		toneFreqParam = new DoubleParameter("Tone Frequency", 1, 192000, "Hz");
		toneFreqParam.setValue(10000d);
		toneFreqParam.addParameterChangeListener(this);
		params.addParameter(toneFreqParam);
		
		samplingFreqParam = new DoubleParameter("Sampling Frequency", 1, 192000, "Hz");
		samplingFreqParam.setValue(48000);
		samplingFreqParam.addParameterChangeListener(this);
		params.addParameter(samplingFreqParam);
		
		fractOffsetParam = new DoubleParameter("Wavelength Fract Offset", 0d, 1d);
		fractOffsetParam.setValue(0d);
		fractOffsetParam.addParameterChangeListener(this);
		params.addParameter(fractOffsetParam);
		
		randOffsetParam = new ButtonParameter("Wavelengh Fract Offset", "Randomize");
		randOffsetParam.addParameterChangeListener(this);
		params.addParameter(randOffsetParam);
		
		GriddedParameterListEditor editor = new GriddedParameterListEditor(params, 1, params.size());
		
		this.add(editor, BorderLayout.NORTH);
		
		gw = new GraphWidget();
		this.add(gw, BorderLayout.CENTER);
		
		updatePlot();
		
		this.setTitle("Sampling GUI");
		this.setVisible(true);
	}
	
	private void updatePlot() {
		double samplingRate = 1d/samplingFreqParam.getValue();
		double toneRate = 1d/toneFreqParam.getValue();
		
		double sampleFractWavelength = samplingRate/toneRate;
		
		double fractOffset = fractOffsetParam.getValue();
		
		System.out.println("sampling rate: "+samplingRate);
		System.out.println("tone rate: "+toneRate);
		System.out.println("sampling fract: "+sampleFractWavelength);
		
		DiscretizedFunc sampledFunc = buildFunc(sampleFractWavelength, fractOffset);
		DiscretizedFunc pureFunc = buildFunc(sampleFractWavelength/1000d, fractOffset);
		
		sampledFunc.setName("Sampled Signal, "+toneFreqParam.getValue().floatValue()
				+" Hz, sampled at "+samplingFreqParam.getValue().floatValue()+" Hz "
						+ "("+(float)(1d/sampleFractWavelength)+" per wavelength)");
		pureFunc.setName("Pure Signal, "+toneFreqParam.getValue().floatValue()+" Hz");
		
		List<PlotElement> funcs = Lists.newArrayList();
		funcs.add(pureFunc);
		funcs.add(sampledFunc);
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GREEN));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, PlotSymbol.FILLED_CIRCLE, 5f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Sampling Test", "Wavelength #", "Y");
		gw.setPlotSpec(spec);
	}
	
	private static DiscretizedFunc buildFunc(double sampleFractWavelength, double fractOffset) {
		if (sampleFractWavelength < 4e-3)
			sampleFractWavelength = 4e-3;
		ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
		
		double angleOffset = fractOffset*Math.PI*2;
		
		for (double x=0; x<=wavelengths_to_plot; x+=sampleFractWavelength) {
			double angle = Math.PI*x*2;
			double y = Math.sin(angle+angleOffset);
			func.set(x, y);
		}
		
		return func;
	}

	@Override
	public void parameterChange(ParameterChangeEvent event) {
		if (event.getParameter() == randOffsetParam) {
			fractOffsetParam.setValue(Math.random());
			fractOffsetParam.getEditor().refreshParamEditor();
		} else {
			updatePlot();
		}
	}
	
	public static void main(String[] args) {
		new SamplingGUI();
		
		/*
		 * can build jar file with this ant in MultiAppBuilder.xml:
		 * 
		 * 
	
	<target name="build.sample.gui">
		<ant antfile="${app.build.file}" target="build.app">
			<property name="app.short.name" value="SamplingGUI" />
			<property name="app.main.class" value="scratch.kevin.SamplingGUI" />
			<property name="javac.includes"     value="scratch/kevin/SamplingGUI.java" />
			<property name="javac.excludes" value="" />
			<property name="resource.target" value="resource.misc.required" />
			<property name="pg" value="true" />
		</ant>
	</target>
		 */
	}

}
