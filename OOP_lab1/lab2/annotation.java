package lab2;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

// https://javarush.com/groups/posts/1896-java-annotacii-chto-ehto-i-kak-ehtim-poljhzovatjhsja
@Retention(RetentionPolicy.RUNTIME) // аннотация которая сохраняется после компиляции и подгружается JVM
@Target({ElementType.METHOD, ElementType.FIELD}) // Тип объекта над которым указывается
public @interface annotation {
    int value();
}
