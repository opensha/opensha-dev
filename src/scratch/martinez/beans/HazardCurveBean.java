package scratch.martinez.beans;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.lang.reflect.InvocationTargetException;
import java.rmi.RemoteException;
import java.util.ArrayList;

import javax.swing.JPanel;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.nshmp.sha.gui.beans.GuiBeanAPI;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.HazardCurveCalculatorAPI;
import org.opensha.sha.earthquake.ERF_Ref;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.gui.beans.ERF_GuiBean;
import org.opensha.sha.gui.beans.IMR_GuiBean;
import org.opensha.sha.gui.beans.IMR_GuiBeanAPI;
import org.opensha.sha.gui.beans.Site_GuiBean;
import org.opensha.sha.imr.ScalarIMR;

public class HazardCurveBean implements GuiBeanAPI, IMR_GuiBeanAPI {
	private JPanel embedVis = null;
	/* The object class names for all the supported Eqk Rup Forecasts */
	private static ArrayList<ERF_Ref> erfs = new ArrayList<ERF_Ref>();
	
	private HazardCurveCalculatorAPI calculator = null;
	private boolean isHazardCalcComplete = true;

	/* A static initialization on startup */
	static {
		erfs.add(ERF_Ref.FRANKEL_ADJUSTABLE_96);
		erfs.add(ERF_Ref.FRANKEL_02);
		erfs.add(ERF_Ref.WGCEP_UCERF_1);
		erfs.add(ERF_Ref.POISSON_FLOATING_FAULT);
		erfs.add(ERF_Ref.POISSON_FAULT);
		erfs.add(ERF_Ref.POINT_SOURCE);
		erfs.add(ERF_Ref.POINT_SOURCE_MULTI_VERT);
		erfs.add(ERF_Ref.WGCEP_02);
	}
	
	////////////////////////////////////////////////////////////////////////////////
	//              Sub-Beans to Drive this Wrapper Convenience Bean              //
	////////////////////////////////////////////////////////////////////////////////
	private Site_GuiBean siteBean = null;
	private ERF_GuiBean erfBean = null;
	private IMR_GuiBean imrBean = null;
	
	////////////////////////////////////////////////////////////////////////////////
	//                             Public Constructors                            //
	////////////////////////////////////////////////////////////////////////////////
	
	/**
	 * Creates a new <code>HazardCurveBean</code> object that can be visualized in any
	 * number of available ways.  The parameters for this bean define the defaults to
	 * use, however these can later be changed.  Any of the parameters can be null, in
	 * such a case, the defaults will be used.
	 */
	public HazardCurveBean(Site_GuiBean siteBean, ERF_GuiBean erfBean, IMR_GuiBean imrBean) {
		// Initialize the site bean
		if(siteBean == null)
			this.siteBean = new Site_GuiBean();
		else
			this.siteBean = siteBean;
		
		// Initialize the erf bean
		if(erfBean == null)
			try {
				this.erfBean = new ERF_GuiBean(erfs);
			} catch (InvocationTargetException ex) {
				System.err.println("Failed creating the erfBean!" + ex.getMessage());
			}
		else
			this.erfBean = erfBean;
		
		// Initialize the imr bean
		if(imrBean == null)
			this.imrBean = new IMR_GuiBean(this);
		else
			this.imrBean = imrBean;
		
		jbInit();
	}
	
	/** A convenience wrapper function to the main constructor */
	public HazardCurveBean() {this(null, null, null);}
	/** A convenience wrapper function to the main constructor */
	public HazardCurveBean(Site_GuiBean siteBean) {this(siteBean, null, null);}
	/** A convenience wrapper function to the main constructor */
	public HazardCurveBean(ERF_GuiBean erfBean) {this(null, erfBean, null);}
	/** A convenience wrapper function to the main constructor */
	public HazardCurveBean(IMR_GuiBean imrBean) {this(null, null, imrBean);}
	/** A convenience wrapper function to the main constructor */
	public HazardCurveBean(Site_GuiBean siteBean, ERF_GuiBean erfBean) {this(siteBean, erfBean, null);}
	/** A convenience wrapper function to the main constructor */
	public HazardCurveBean(Site_GuiBean siteBean, IMR_GuiBean imrBean) {this(siteBean, null, imrBean);}
	/** A convenience wrapper function to the main constructor */
	public HazardCurveBean(ERF_GuiBean erfBean, IMR_GuiBean imrBean) {this(null, erfBean, imrBean);}
	
	////////////////////////////////////////////////////////////////////////////////
	//                   Public Functions to Help Use the Bean                    //
	////////////////////////////////////////////////////////////////////////////////
	
	/* Public Utility Fucntions */
	
	public DiscretizedFunc computeHazardCurve(ArrayList<Double> imls) {
		isHazardCalcComplete = false;
		DiscretizedFunc newFunc = new ArbitrarilyDiscretizedFunc();
		for(int i = 0; i < imls.size(); ++i)
			newFunc.set(imls.get(i), 0.0);
		
		try {
			newFunc = calculator.getHazardCurve(newFunc, siteBean.getSite(), 
					imrBean.getSelectedIMR_Instance(), (ERF) erfBean.getSelectedERF());
		} catch (InvocationTargetException e) {
			e.printStackTrace();
		}
		

		isHazardCalcComplete = true;
		return newFunc;
	}
	
	/**
	 * Allows an observer to ask whether or not its application is complete yet.
	 * 
	 * @return True if the hazard curve computations have finished, false otherwise
	 */
	public boolean isCalculationFinished() {
		return isHazardCalcComplete;
	}
	
	
	/* Public Getters */
	
	/**
	 * Get the underlying <code>Site_GuiBean</code>.
	 */
	public Site_GuiBean getSiteBean() {return siteBean;}
	/**
	 * Get the underlying <code>IMR_GuiBean</code>.
	 */
	public IMR_GuiBean getImrBean() {return imrBean;}
	/**
	 * Get the underlying <code>ERF_GuiBean</code>.
	 */
	public ERF_GuiBean getErfBean() {return erfBean;}
	
	/* Public Setters */
	
	/**
	 * Set the underlying <code>Site_GuiBean</code>
	 * @param newBean The new bean to use.
	 */
	public void setSiteBean(Site_GuiBean newBean) {siteBean = newBean;}
	/**
	 * Set the underlying <code>IMR_GuiBean</code>
	 * @param newBean The new bean to use.
	 */
	public void setImrBean(IMR_GuiBean newBean) {imrBean = newBean;}
	/**
	 * Set the underlying <code>ERF_GuiBean</code>
	 * @param newBean The new bean to use.
	 */
	public void setErfBean(ERF_GuiBean newBean) {erfBean = newBean;}
	
	////////////////////////////////////////////////////////////////////////////////
	//                 Functions to Implement Interfaces Specified                //
	////////////////////////////////////////////////////////////////////////////////
	/** See the general contract in GuiBeanAPI */
	public Object getVisualization(int type) throws IllegalArgumentException {
		if(type == EMBED)
			return getEmbeddableVisualization();
		else
			throw new IllegalArgumentException("The given type is not currently supported.");
	}

	/** See the general contract in GuiBeanAPI */
	public String getVisualizationClassName(int type) {
		if(type == EMBED)
			return "javax.swing.JPanel";
		else
			return null;
	}

	/** See the general contract in GuiBeanAPI */
	public boolean isVisualizationSupported(int type) {
		return (type == EMBED);
	}
	
	/** See the general contract in IMR_GuiBeanAPI */
	public void updateIM() {
		// Do nothing for now...is this not important?
	}

	/** See the general contract in IMR_GuiBeanAPI */
	public void updateSiteParams() {
		ScalarIMR imr = imrBean.getSelectedIMR_Instance();
		siteBean.replaceSiteParams(imr.getSiteParamsIterator());
		siteBean.validate();
		siteBean.repaint();
	}
	
	////////////////////////////////////////////////////////////////////////////////
	//                      Private Functions Used Internally                     //
	////////////////////////////////////////////////////////////////////////////////

	/* Initializes the bean */
	public void jbInit() {
		calculator = new HazardCurveCalculator();
	}
	
	/* Creates the embeddable JFrame to an EMBED visualization */
	private JPanel getEmbeddableVisualization() {
		if(embedVis == null) {
			embedVis = new JPanel(new GridBagLayout());
			embedVis.add(imrBean, new GridBagConstraints(0, 0, 1, 2, 1.0, 1.0,
					GridBagConstraints.CENTER, GridBagConstraints.VERTICAL,
					new Insets(2, 2, 2, 2), 2, 2));
			embedVis.add(erfBean, new GridBagConstraints(1, 0, 1, 3, 1.0, 1.0,
					GridBagConstraints.CENTER, GridBagConstraints.VERTICAL,
					new Insets(2, 2, 2, 2), 2, 2));
			embedVis.add(siteBean, new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0,
					GridBagConstraints.CENTER, GridBagConstraints.VERTICAL,
					new Insets(2, 2, 2, 2), 2, 2));
		}
		return embedVis;
	}
}
