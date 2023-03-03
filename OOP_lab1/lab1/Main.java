import common.Point;
import strateg.*;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;

public class Main {
    public static void main(String[] args) {
        HashMap<String, MoveStrategy> choosingInStream = new HashMap<>();
        choosingInStream.put("fly", new Fly());
        choosingInStream.put("ride", new RideAHorse());
        choosingInStream.put("stay", new heStands());

        System.out.println("Hello, Write point of destination in format \"x.x x.x\" ");
        Hero hero = new Hero(new Point(0, 0)); // default position

        boolean exit = false;                        // flag for exit game
        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        String[] list;
        while (!exit) {
            hero.move();
            System.out.printf("Now your hero in %f %f %n", hero.getPos().x, hero.getPos().y);

            try {
                list = reader.readLine().split(" "); /* dividing the array into
                                                            elements separated by a space */
                }
            catch (IOException e) {
                list = new String[]{""};   // filling in the null element of the array
            }

            try {
                switch (list[0]) {
                    case "exit" -> exit = true;
                    case "dest" -> hero.setDest(new Point(Float.parseFloat(list[1]), Float.parseFloat(list[2])));
                    case "strategy" -> hero.setMoveStrategy(choosingInStream.get(list[1])); /* selecting a class to pass
                                                                                             an argument to the parent constructor */
                }
            }
            catch(NumberFormatException | IndexOutOfBoundsException e){
                System.out.println("you are invalid");
            }
        }
    }
}