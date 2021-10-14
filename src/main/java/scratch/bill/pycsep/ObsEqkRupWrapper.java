package scratch.bill.pycsep;


import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.time.LocalDateTime;
import java.util.GregorianCalendar;
import java.util.TimeZone;

import org.opensha.commons.geo.Location;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;

public class ObsEqkRupWrapper {

    public static ObsEqkRupList readCsepAscii(String fileName, boolean hasHeader) throws IOException {
        BufferedReader inputReader = new BufferedReader(new FileReader(fileName));
        ObsEqkRupList catalog = new ObsEqkRupList();
        String line;
        String[] split;
        boolean headerLine = true;
        // read and parse files
        while(inputReader.ready()) {
            // split on ',' and remove whitespace
            line = inputReader.readLine();
            line = line.trim();
            split = line.split(",");
            // handle header line
            if (hasHeader && headerLine) {
                headerLine = false;
                continue;
            }
            // first parse fields from csv file
            double longitude = Double.parseDouble(split[0]);
            double latitude = Double.parseDouble(split[1]);
            double magnitude = Double.parseDouble(split[2]);
            LocalDateTime ldt = LocalDateTime.parse(split[3]);
            double depth = Double.parseDouble(split[4]);
            String eventId = split[6];
            // now, convert them to java objects
            GregorianCalendar cal = new GregorianCalendar(TimeZone.getTimeZone("GMT-0:00"));
            cal.clear();
            cal.set(ldt.getYear(), ldt.getMonthValue()-1, ldt.getDayOfMonth(), ldt.getHour(), ldt.getMinute(), ldt.getSecond());
            Location hypoLoc = new Location(latitude, longitude, depth);
            // add rupture to rupturelist
            catalog.add(new ObsEqkRupture(eventId, cal.getTimeInMillis(), hypoLoc, magnitude));

        }
        return catalog;
    }

    public static void writeCsepAscii(String fileName, boolean writeHeader) {
        // Todo: implement catalog writer in pyCSEP format
    }
}
