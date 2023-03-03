import common.Point;
import strateg.MoveStrategy;
import strateg.heStands;


public class Hero {
    private Point pos;
    private Point dest;
    private MoveStrategy moveStrategy;

    public Hero(Point pos){
        this.pos = pos;
        this.dest = pos;
        this.moveStrategy = new heStands();
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

    public void setMoveStrategy(MoveStrategy strategy) {
        this.moveStrategy = strategy;
    }

    public void move() {
        setPos(moveStrategy.getNewPos(pos, dest));
    }
}
