package scratch.kevin.cybershake;


public class CBSiteAmpCalc {

//	public static void main(String[] args) throws DocumentException, InvocationTargetException, IOException {
////		boolean ddwCorr = false;
////		File outputDir = new File("/home/kevin/CyberShake/MCER/gmpe_site_amp_noddw");
////		boolean ddwCorr = true;
////		File outputDir = new File("/home/kevin/CyberShake/MCER/gmpe_site_amp");
////		boolean ddwCorr = false;
////		File outputDir = new File("/home/kevin/CyberShake/MCER/gmpe_site_amp_bssa");
////		boolean ddwCorr = false;
////		File outputDir = new File("/home/kevin/CyberShake/MCER/gmpe_site_amp_all_classes_bssa");
////		File outputDir = new File("/home/kevin/CyberShake/MCER/gmpe_site_amp_all_classes_bssa_redo");
//		boolean ddwCorr = false;
////		File outputDir = new File("/home/kevin/CyberShake/MCER/gmpe_site_amp_all_classes_mean_redo");
//		boolean twoPercentIn50 = false;
//		boolean noDeterministic = true;
//		File outputDir = new File("/home/kevin/CyberShake/MCER/"
//				+ "gmpe_site_amp_all_classes_mean_no_determ");
////		File outputDir = new File("/tmp/asdf");
//		
//		if (twoPercentIn50)
//			outputDir = new File(outputDir.getParentFile(), outputDir.getName()+"_2pin50");
//		
//		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
//		
////		AttenuationRelationship meanGMPE = AttenRelRef.BSSA_2014.instance(null);
////		AttenuationRelationship meanGMPE = new NGAWest_2014_Averaged_AttenRel(null, false);
//		
//		List<Integer> runIDs = Lists.newArrayList(
////				2657, 3037, 2722, 3022, 3030, 3027, 2636,
////				2638, 2660, 2703, 3504, 2988, 2965, 3007);
//				
//				2988, 2657);
//				
////				2657, 3037, 2722, 3022, 3030, 3027, 2636);
////				2638, 2660, 2703, 3504, 2988, 2965, 3007);
//		
////				2636);
//		
//		File asceFile = new File("/home/kevin/CyberShake/MCER/ASCE7-10_Sms_Sm1_TL_det LL for 14 sites.xls");
//		
//		boolean backSeis = true;
//		
//		ERF probERF = MeanUCERF2_ToDB.createUCERF2ERF();
//		ERF detERF;
//		if (backSeis) {
//			probERF.setParameter(UCERF2.BACK_SEIS_NAME, UCERF2.BACK_SEIS_INCLUDE);
//			detERF = MeanUCERF2_ToDB.createUCERF2ERF();
//			probERF.setParameter(MeanUCERF2.CYBERSHAKE_DDW_CORR_PARAM_NAME, ddwCorr);
//			detERF.setParameter(MeanUCERF2.CYBERSHAKE_DDW_CORR_PARAM_NAME, ddwCorr);
//			probERF.updateForecast();
//			detERF.updateForecast();
//		} else {
//			probERF.setParameter(MeanUCERF2.CYBERSHAKE_DDW_CORR_PARAM_NAME, ddwCorr);
//			probERF.updateForecast();
//			detERF = probERF;
//		}
//		List<AttenuationRelationship> attenRels = Lists.newArrayList();
//
//		String attenFiles = "src/org/opensha/sha/cybershake/conf/cb2014.xml,src/org/opensha/sha/cybershake/conf/cy2014.xml"
//				+ ",src/org/opensha/sha/cybershake/conf/bssa2014.xml,src/org/opensha/sha/cybershake/conf/ask2014.xml";
//		for (String attenRelFile : HazardCurvePlotter.commaSplit(attenFiles)) {
//			AttenuationRelationship attenRel = AttenRelSaver.LOAD_ATTEN_REL_FROM_FILE(attenRelFile);
//			attenRels.add(attenRel);
//		}
//		Preconditions.checkArgument(!attenRels.isEmpty(), "Must specify at least 1 GMPE");
//		MultiIMR_Averaged_AttenRel meanGMPE = new MultiIMR_Averaged_AttenRel(attenRels);
//		
//		meanGMPE.setParamDefaults();
//		List<AttenuationRelationship> meanGMPEList = Lists.newArrayList();
//		meanGMPEList.add(meanGMPE);
//
//		List<Double> periods = Lists.newArrayList(0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,
//				0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,7.5,10.0);
////		List<Double> periods = Lists.newArrayList(2.0,3.0,4.0,5.0,7.5,10.0);
//		
////		double[] vs30s = {1620d};
////		boolean[] nullBasins = {true};
//		double[] vs30s = {1620d, 1524d, 914d, 762d, 488d, 366d, 265d, 183d, 155d};
//		boolean[] nullBasins = {true, false};
////		double[] vs30s = {1524d, 762d, 488d, 265d};
////		boolean[] nullBasins = {true};
//		
//		DBAccess db = Cybershake_OpenSHA_DBApplication.getDB();
//		Runs2DB runs2db = new Runs2DB(db);
//		CybershakeSiteInfo2DB sites2db = new CybershakeSiteInfo2DB(db);
//		CyberShakeComponent comp = CyberShakeComponent.RotD100;
//		double percentile = GMPEDeterministicComparisonCalc.default_percentile;
//		
////		HSSFSheet asceSheet = null;
////		FormulaEvaluator evaluator = null;
////		if (includeASCE) {
////			HSSFWorkbook wb;
////			try {
////				POIFSFileSystem fs = new POIFSFileSystem(new FileInputStream(asceFile));
////				wb = new HSSFWorkbook(fs);
////			} catch (Exception e1) {
////				System.err.println("Couldn't load input file. Make sure it's an xls file and NOT an xlsx file.");
////				throw ExceptionUtils.asRuntimeException(e1);
////			}
////			asceSheet = wb.getSheetAt(0);
////			evaluator = wb.getCreationHelper().createFormulaEvaluator();
////		}
//		
//		List<CybershakeIM> forceAddIMs = Lists.newArrayList();
//		for (Double period : periods)
//			if (period < 2d || period == 6d)
//				forceAddIMs.add(new CybershakeIM(-1, IMType.SA, period, null, comp));
//		
//		for (int runID : runIDs) {
//			System.out.println("**************** RUN ID: "+runID);
//			
//			CybershakeRun run = runs2db.getRun(runID);
//			CybershakeSite site = sites2db.getSiteFromDB(run.getSiteID());
//			
//			GMPEDeterministicComparisonCalc detCalc = new GMPEDeterministicComparisonCalc(run, site, comp, periods,
//					percentile, detERF, attenRels, null);
//			
//			int velModelID = run.getVelModelID();
//			OrderedSiteDataProviderList providers = HazardCurvePlotter.createProviders(velModelID);
//			List<SiteDataValue<?>> origSiteDatas = providers.getBestAvailableData(site.createLocation());
//			
//			List<List<SiteDataValue<?>>> siteDatasList = Lists.newArrayList();
//			
//			// first scenario, normal site data
//			siteDatasList.add(origSiteDatas);
////			// null basin, Wills Vs30
////			siteDatasList.add(getReplaced(origSiteDatas,
////					new SiteDataValue<Double>(SiteData.TYPE_DEPTH_TO_1_0, SiteData.TYPE_FLAG_INFERRED, Double.NaN),
////					new SiteDataValue<Double>(SiteData.TYPE_DEPTH_TO_2_5, SiteData.TYPE_FLAG_INFERRED, Double.NaN)));
////			// null basin, rock
////			siteDatasList.add(getReplaced(origSiteDatas,
////					new SiteDataValue<Double>(SiteData.TYPE_DEPTH_TO_1_0, SiteData.TYPE_FLAG_INFERRED, Double.NaN),
////					new SiteDataValue<Double>(SiteData.TYPE_DEPTH_TO_2_5, SiteData.TYPE_FLAG_INFERRED, Double.NaN),
////					new SiteDataValue<Double>(SiteData.TYPE_VS30, SiteData.TYPE_FLAG_INFERRED, 760d)));
////			// CVM basin, rock
////			siteDatasList.add(getReplaced(origSiteDatas,
////					new SiteDataValue<Double>(SiteData.TYPE_VS30, SiteData.TYPE_FLAG_INFERRED, 760d)));
////			// null basin, soft
////			siteDatasList.add(getReplaced(origSiteDatas,
////					new SiteDataValue<Double>(SiteData.TYPE_DEPTH_TO_1_0, SiteData.TYPE_FLAG_INFERRED, Double.NaN),
////					new SiteDataValue<Double>(SiteData.TYPE_DEPTH_TO_2_5, SiteData.TYPE_FLAG_INFERRED, Double.NaN),
////					new SiteDataValue<Double>(SiteData.TYPE_VS30, SiteData.TYPE_FLAG_INFERRED, 155d)));
////			// CVM basin, soft
////			siteDatasList.add(getReplaced(origSiteDatas,
////					new SiteDataValue<Double>(SiteData.TYPE_VS30, SiteData.TYPE_FLAG_INFERRED, 155d)));
//			
//			for (boolean nullBasin : nullBasins) {
//				for (double vs30 : vs30s) {
//					if (nullBasin)
//						siteDatasList.add(getReplaced(origSiteDatas,
//								new SiteDataValue<Double>(SiteData.TYPE_DEPTH_TO_1_0, SiteData.TYPE_FLAG_INFERRED, Double.NaN),
//								new SiteDataValue<Double>(SiteData.TYPE_DEPTH_TO_2_5, SiteData.TYPE_FLAG_INFERRED, Double.NaN),
//								new SiteDataValue<Double>(SiteData.TYPE_VS30, SiteData.TYPE_FLAG_INFERRED, vs30)));
//					else
//						siteDatasList.add(getReplaced(origSiteDatas,
//								new SiteDataValue<Double>(SiteData.TYPE_VS30, SiteData.TYPE_FLAG_INFERRED, vs30)));
//				}
//			}
//			
//			CSVFile<String> csv = new CSVFile<String>(true);
//			
//			List<String> vs30Line = Lists.newArrayList();
//			List<String> z10Line = Lists.newArrayList();
//			List<String> z25Line = Lists.newArrayList();
//			vs30Line.add("Vs30");
//			z10Line.add("Z1.0");
//			z25Line.add("Z2.5");
//			List<List<String>> periodLines = Lists.newArrayList();
//			for (List<SiteDataValue<?>> siteDatas : siteDatasList) {
//				vs30Line.add(getDataStr(siteDatas, SiteData.TYPE_VS30));
//				z10Line.add(getDataStr(siteDatas, SiteData.TYPE_DEPTH_TO_1_0));
//				z25Line.add(getDataStr(siteDatas, SiteData.TYPE_DEPTH_TO_2_5));
//			}
//			csv.addLine(vs30Line);
//			csv.addLine(z10Line);
//			csv.addLine(z25Line);
//			
//			for (Double period : periods) {
//				List<String> line = Lists.newArrayList();
//				line.add(period+"");
//				for (int i=0; i<siteDatasList.size(); i++)
//					line.add("");
//				periodLines.add(line);
//				csv.addLine(line);
//			}
//			
//			for (int s=0; s<siteDatasList.size(); s++) {
//				List<SiteDataValue<?>> siteDatas = siteDatasList.get(s);
//				for (AttenuationRelationship attenRel : attenRels)
//					attenRel.setParamDefaults();
//				
//				RTGMCalc rtgmCalc = new RTGMCalc(runID, comp, null, db);
//				rtgmCalc.setUse2PercentIn50(twoPercentIn50);
//				rtgmCalc.setSiteDatas(siteDatas);
//				rtgmCalc.setGMPEs(probERF, meanGMPEList);
//				rtgmCalc.setForceAddIMs(forceAddIMs);
//				Preconditions.checkState(rtgmCalc.calc());
//				
//				DiscretizedFunc probFunc = MCErCalcUtils.saToPsuedoVel(rtgmCalc.getGMPESpectrumMap().get(comp).get(0));
//				
//				DiscretizedFunc detFunc = null;
//				if (!twoPercentIn50 && !noDeterministic) {
//					detCalc.setSiteData(siteDatas);
//					
//					detCalc.calc();
//					
//					Table<Double, AttenuationRelationship, DeterministicResult> detVals = detCalc.getResults();
//					detFunc = new ArbitrarilyDiscretizedFunc();
//					for (double period : periods) {
//						double maxY = 0;
//						for (AttenuationRelationship attenRel : attenRels)
//							maxY = Math.max(maxY, detVals.get(period, attenRel).getVal());
//						detFunc.set(period, maxY);
//					}
//					detFunc = MCErCalcUtils.saToPsuedoVel(detFunc);
//					
//					Preconditions.checkState(probFunc.size() == detFunc.size(), probFunc.size()+" != "+detFunc.size());
//				}
//				
//				double vs30 = (Double)attenRels.get(0).getParameter(Vs30_Param.NAME).getValue();
//				
//				DiscretizedFunc asceDeterm = MCErCalcUtils.saToPsuedoVel(
//						MCERDataProductsCalc.calc(probFunc.deepClone(), vs30, site.createLocation()));
////				if (includeASCE) {
////					DiscretizedFunc xVals = detFunc.deepClone();
////					HSSFRow row = null;
////					for (int r=0; r<=asceSheet.getLastRowNum(); r++) {
////						HSSFRow testRow = asceSheet.getRow(r);
////						HSSFCell nameCell = testRow.getCell(0);
////						if (nameCell != null && nameCell.getStringCellValue().trim().equals(site.short_name)) {
////							row = testRow;
////							break;
////						}
////					}
////					Preconditions.checkState(row != null, "Couldn't find site "+site.short_name+" in ASCE spreadsheet");
////					double tl = MCERDataProductsCalc.loadASCEValue(row.getCell(4), evaluator);
////					double detASCE =  MCERDataProductsCalc.loadASCEValue(row.getCell(7), evaluator);
////					asceDeterm = MCERDataProductsCalc.calcASCE(xVals, detASCE, tl);
////				}
//				
//				DiscretizedFunc mcer = MCERDataProductsCalc.calcMCER(detFunc, probFunc, asceDeterm);
//				Preconditions.checkState(mcer.size() == periods.size());
//				
//				for (int i=0; i<periods.size(); i++) {
//					double period = periods.get(i);
//					double val = mcer.getY(period);
//					periodLines.get(i).set(s+1, val+"");
//				}
//			}
//			String name = site.short_name+"_run"+run.getRunID()+"_GMPE_SiteAmp_"+comp.getShortName()+".csv";
//			
//			File outputFile = new File(outputDir, name);
//			System.out.println("**************** WRITING RESULTS TO "+outputFile.getAbsolutePath());
//			csv.writeToFile(outputFile);
//			
//			// now write SA version
//			CSVFile<String> saCSV = new CSVFile<String>(true);
//			// headers
//			saCSV.addLine(csv.getLine(0));
//			saCSV.addLine(csv.getLine(1));
//			saCSV.addLine(csv.getLine(2));
//			for (int i=saCSV.getNumRows(); i<csv.getNumRows(); i++) {
//				List<String> line = Lists.newArrayList(csv.getLine(i));
//				double period = Double.parseDouble(line.get(0));
//				for (int j=1; j<line.size(); j++) {
//					double val = Double.parseDouble(line.get(j));
//					val *= 2*Math.PI/(period*HazardCurveComputation.CONVERSION_TO_G);
//					line.set(j, val+"");
//				}
//				saCSV.addLine(line);
//			}
//			
//			name = site.short_name+"_run"+run.getRunID()+"_GMPE_SiteAmp_"+comp.getShortName()+"_sa.csv";
//			
//			outputFile = new File(outputDir, name);
//			System.out.println("**************** WRITING RESULTS TO "+outputFile.getAbsolutePath());
//			saCSV.writeToFile(outputFile);
//		}
//		
//		System.exit(0);
//	}
//	
//	private static List<SiteDataValue<?>> getReplaced(List<SiteDataValue<?>> origSiteDatas,
//			SiteDataValue<Double>... replacements) {
//		List<SiteDataValue<?>> replaced = Lists.newArrayList();
//		
//		int numReplaced = 0;
//		
//		valLoop:
//		for (SiteDataValue<?> orig : origSiteDatas) {
//			SiteDataValue<?> val = orig;
//			for (SiteDataValue<Double> replacement : replacements) {
//				if (replacement.getDataType().equals(val.getDataType())) {
//					val = replacement;
//					numReplaced++;
//					if (replacement.getValue().isNaN())
//						// don't add this one, we want null
//						continue valLoop;
//					break;
//				}
//			}
//			replaced.add(val);
//		}
//		Preconditions.checkState(numReplaced == replacements.length);
//		
//		return replaced;
//	}
//	
//	private static SiteDataValue<?> getByType(List<SiteDataValue<?>> vals, String type) {
//		for (SiteDataValue<?> val : vals)
//			if (val.getDataType().equals(type))
//				return val;
//		return null;
//	}
//	
//	public static String getDataStr(List<SiteDataValue<?>> vals, String type) {
//		SiteDataValue<?> val = getByType(vals, type);
//		if (val != null)
//			return val.getValue().toString();
//		return "(null)";
//	}

}
