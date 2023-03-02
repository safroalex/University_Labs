import startHero.Point;
import strateg.heStands;
import strateg.start_stop_point;

public class Hero {
    private Point pos;
    private Point dest;
    private start_stop_point startStopPoint;

    public Hero(Point pos){
        this.pos = pos;
        this.dest = pos;
        this.startStopPoint = new heStands();
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

    public void setMoveStrategy(start_stop_point startStopPoint) {
        this.startStopPoint = startStopPoint;
    }

    public void move() {
        setPos(startStopPoint.newPos(pos, dest));
    }
}
