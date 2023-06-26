const int RELAY_ENABLE = 5;

void setup() {
  // put your setup code here, to run once:
  pinMode(RELAY_ENABLE,OUTPUT);
  Serial.begin(115200);

}

void loop() {
  // put your main code here, to run repeatedly:
  Serial.println("Relay On");
  digitalWrite(RELAY_ENABLE, LOW); 
  delay(1000);

  Serial.println("Relay OFF");
  digitalWrite(RELAY_ENABLE,HIGH);
  delay(1000);

}
