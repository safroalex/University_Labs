package lab2;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

public class Main {
    public static void main(String[] args) {
        classForCalling classForCalling = new classForCalling();
        try {
            callOurMethods(classForCalling);
        }
        catch (InvocationTargetException | IllegalAccessException e) {
            e.printStackTrace();
        }
    }

    private static void callOurMethods(Object object) throws InvocationTargetException, IllegalAccessException {
        Class clazz = object.getClass();
        for (Method method : clazz.getDeclaredMethods()) {
            if (method.isAnnotationPresent(annotation.class)) {
                for (int i = 0; i < method.getAnnotation(annotation.class).value(); i++) {
                    method.setAccessible(true);
                    method.invoke(object);
                }
            }
        }

    }
}