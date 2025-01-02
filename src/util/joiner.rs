use std::fmt::{Display, Formatter};

pub(crate) struct Joiner<'a, 'b, T: Display> {
    sep: &'a str,
    list: &'b [T],
}

impl<'a, 'b, T: Display> Joiner<'a, 'b, T> {
    pub(crate) fn new(sep: &'a str, list: &'b [T]) -> Joiner<'a, 'b, T> {
        Joiner { sep, list }
    }
}

impl<T: Display> Display for Joiner<'_, '_, T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let iter = self.list.iter();
        write_iter_fmt(f, iter, self.sep)?;
        Ok(())
    }
}

fn join_iter<T, I, E, P>(sep: &str, mut iter: I, mut print_item: P) -> Result<(), E>
where I: Iterator<Item=T>, P: FnMut(&str, T) -> Result<(), E>
{
    if let Some(first) = iter.next() {
        print_item("", first)?;
        for item in iter {
            print_item(sep, item)?;
        }
    }
    Ok(())
}

pub(crate) fn write_iter_fmt<T, I, W>(w: &mut W, iter: I, sep: &str) -> std::fmt::Result
where T: Display, I: Iterator<Item=T>, W: std::fmt::Write
{
    join_iter(sep, iter, |sep, item| write!(w, "{}{}", sep, item))
}

pub(crate) fn write_iter_io<T, I, W>(w: &mut W, iter: I, sep: &str) -> std::io::Result<()>
where T: Display, I: Iterator<Item=T>, W: std::io::Write
{
    join_iter(sep, iter, |sep, item| write!(w, "{}{}", sep, item))
}
