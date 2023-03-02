package lab1;

import lab1.common.Point;
import lab1.strategies.MoveStrategy;
import lab1.strategies.Stay;

public class Hero {
    private Point pos;
    private Point dest;
    private MoveStrategy moveStrategy;

    public Hero(Point pos) {
        this.pos = pos;
        this.dest = pos;
        this.moveStrategy = new Stay();
    }

    public void setPos(Point pos) {
        this.pos = pos;
    }

    public Point getPos() {
        return pos;
    }

    public void setDest(Point destination) {
        this.dest = destination;
    }

    public void setMoveStrategy(MoveStrategy moveStrategy) {
        this.moveStrategy = moveStrategy;
    }

    public void move() {
        setPos(moveStrategy.getNewPos(pos, dest));
    }
}
