package scratch.aftershockStatistics.cmu;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

/**
 * Created by clark on 12/12/2016.
 */
public class Configuration {

    private String db_host = "127.0.0.1";
    private int db_port = 27017;
    private String db_name = "usgs";
    private String db_user = "usgs";
    private String db_password = "usgs";
    private String activemq_host = "127.0.0.1";
    private String activemq_port = "61616";
    private String activemq_user = "admin";
    private String activemq_password = "admin";
    private String path = "config.properties";

    public Configuration() {
        if (System.getProperty("path") != null) {
            path = System.getProperty("path");
        }

        Properties prop = new Properties();

        //InputStream inputStream = getClass().getClassLoader().getResourceAsStream(path);
        InputStream inputStream = null;
        try {
            inputStream = new FileInputStream(path);
        } catch (FileNotFoundException e) {
            System.out.println("File not found at " + path);
            System.exit(0);
        }

        if(inputStream != null){
            try {
                prop.load(inputStream);
            } catch (IOException e) {
                System.out.println("Error while reading file at " + path);
                System.exit(0);
            }
        } else {
            System.out.println("File at " + path + " not found or it can't be read");
            System.exit(0);
        }

        try {
            this.db_host = prop.getProperty("db_host");
            this.db_port = Integer.parseInt(prop.getProperty("db_port"));
            this.db_name = prop.getProperty("db_name");
            this.db_user = prop.getProperty("db_user");
            this.db_password = prop.getProperty("db_password");
            this.activemq_host = prop.getProperty("activemq_host");
            this.activemq_port = prop.getProperty("activemq_port");
            this.activemq_user = prop.getProperty("activemq_user");
            this.activemq_password = prop.getProperty("activemq_password");
        } catch(Exception e) {
            System.out.println("Error reading properties from file: " + path);
            System.exit(0);
        }
    }

    public String getDb_password() {
        return db_password;
    }

    public String getDb_host() {
        return db_host;
    }

    public int getDb_port() {
        return db_port;
    }

    public String getDb_name() {
        return db_name;
    }

    public String getDb_user() {
        return db_user;
    }

    public String getActivemq_host() {
        return activemq_host;
    }

    public String getActivemq_port() {
        return activemq_port;
    }

    public String getActivemq_user() {
        return activemq_user;
    }

    public String getActivemq_password() {
        return activemq_password;
    }
}
