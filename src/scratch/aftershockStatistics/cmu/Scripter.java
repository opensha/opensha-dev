package scratch.aftershockStatistics.cmu;

import javax.jms.*;

import org.apache.activemq.ActiveMQConnectionFactory;
import javax.jms.MessageConsumer;

public final class Scripter {

  public static void main(final String[] args) throws Exception {
      Configuration config = new Configuration();

  	System.out.println("Scripter Service Starting...");
    final ConnectionFactory connectionFactory = new ActiveMQConnectionFactory("tcp://" + config.getActivemq_host() + ":" + config.getActivemq_port());
    final Connection connection = connectionFactory.createConnection(config.getActivemq_user(), config.getActivemq_password());
    
    Runtime.getRuntime().addShutdownHook(new Thread(){
        @Override
        public void run(){
            try {
            	System.out.println("Shutdown Hook Activated...");
                if(connection != null){
                    connection.close();
                }
            } catch (JMSException e) {
                e.printStackTrace();
            }
        }
    });
    
    try{
      final Session session = connection.createSession(false, Session.AUTO_ACKNOWLEDGE);

        /** PDL Consumer */
        final Queue PDLQueue = session.createQueue("PDL");
      MessageConsumer PDLConsumer = session.createConsumer(PDLQueue);
        PDLConsumer.setMessageListener(new ComcatListener(session));

        /** Comcat Consumer */
      final Queue ComcatQueue = session.createQueue("Comcat");
      MessageConsumer ComcatConsumer = session.createConsumer(ComcatQueue);
        ComcatConsumer.setMessageListener(new ComcatListener(session));
        connection.start();

        System.out.println("Scripter Service Listening to the Queues PDL and Comcat...");
      while(true){
          Thread.sleep(Integer.MAX_VALUE);
      }
    } finally {
    	System.out.println("Scripter Service Exiting...");
      if (connection != null) {
          connection.close();
      }
    }
  }
}