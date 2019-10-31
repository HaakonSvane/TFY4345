

def foo(f_external = lambda x: 0):
    return f_external(10)

print(foo())