use std::fmt::Display;

pub(crate) struct DisplayOption<T: Display> {
    option: Option<T>
}

impl<T: Display> DisplayOption<T> {
    pub(crate) fn new(option: Option<T>) -> DisplayOption<T> {
        DisplayOption { option }
    }
}

impl<T: Display> Display for DisplayOption<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match &self.option {
            Some(value) => write!(f, "{}", value),
            None => write!(f, "None"),
        }
    }
}