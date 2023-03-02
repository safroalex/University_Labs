package lab1;

import lab1.common.Point;
import lab1.strategies.*;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;

public class Main {
    public static void main(String[] args) {
        HashMap<String, MoveStrategy> moveStrategyMap = new HashMap<>();
        moveStrategyMap.put("stay", new Stay());
        moveStrategyMap.put("fly", new Fly());
        moveStrategyMap.put("ride", new RideAHorse());
        moveStrategyMap.put("walk", new Walk());
        moveStrategyMap.put("teleport", new Teleport());

        Hero hero = new Hero(new Point(0, 0));
        boolean exitGameLoop = false;
        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        String[] args_ = {};
        while (!exitGameLoop) {
            hero.move();
            System.out.printf("Hero position is (%f, %f)%n", hero.getPos().x, hero.getPos().y);
            try {
                args_ = reader.readLine().split(" ");
            } catch (IOException e) {
                args_ = new String[]{""};
            }
            if (args_.length == 0) {
                args_ = new String[]{"pass"};
            }
            try {
                switch (args_[0]) {
                    case "exit": {
                        exitGameLoop = true;
                        break;
                    }
                    case "dest": {
                        hero.setDest(new Point(Float.parseFloat(args_[1]), Float.parseFloat(args_[2])));
                        break;
                    }
                    case "strategy": {
                        hero.setMoveStrategy(moveStrategyMap.get(args_[1]));
                        break;
                    }
                }
            } catch (NumberFormatException | IndexOutOfBoundsException e) {
                System.out.println("invalid args for command");
            }
        }
    }
}