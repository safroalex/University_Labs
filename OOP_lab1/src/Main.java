import startHero.Point;
import strateg.Fly;
import strateg.RideAHorse;
import strateg.heStands;
import strateg.start_stop_point;

import java.util.HashMap;

public class Main {
    public static void main(String[] args) {
        HashMap<String, start_stop_point> choosingInStream = new HashMap<>();
        choosingInStream.put("fly", new Fly());
        choosingInStream.put("ride", new RideAHorse());
        choosingInStream.put("stay", new heStands());

        Hero hero = new Hero(new Point(0, 0)); // позиция по умолчанию
    }
}